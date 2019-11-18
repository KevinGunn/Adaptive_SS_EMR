set search_path to mimiciii;

/*
This code extracts all variables used in this model.
This includes pre HE serum creatinine, post HE serum creatinine, mean MAP and total urine output for 3 hours previous to onset of HE,
obtain age, gender, vital signs and service type (ICU), elixhauser scores, avg. temperature in 3 hour window before HE, and SAPS I scores.
Use table "hypo_cohort_fluids_vp" as starting point. This table includes information about hypotensive episode, treatment, database source, and length of stay (los).
*/

select count(*) from hypo_cohort_fluids_vp_subjects;

select * from d_labitems where lower(label) like '%creatinine%';
--5,705 patients in data set.

-- Step 1. Get Serum Creatinine measures

/* Get baseline serum creatinine */
with pre24creat AS (
  SELECT h.*, l.charttime as pre_creat_time, l.valuenum as pre_creat, (h.he_onset - l.charttime) as first_creat_length
    FROM hypo_cohort_fluids_vp_subjects h
    LEFT JOIN labevents l ON
    l.SUBJECT_ID = h.SUBJECT_ID 
    WHERE l.itemid = 50912 AND  ( h.he_onset - l.charttime ) between '00:00:00' and '24:00:00'
    order by SUBJECT_ID
)
, last_creat AS(
	select distinct subject_id, icustay_id,
    	   last_value( pre_creat_time ) over (partition by subject_id, icustay_id) as last_creat_time
    from pre24creat
)
select p.* 
into hypo_cohort_basecreat
from pre24creat p, last_creat l
where p.icustay_id=l.icustay_id and p.pre_creat_time=l.last_creat_time;

/* Get post-episode serum creatinine */
with post72creat AS (
  SELECT h.*, l.charttime as post_creat_time, l.valuenum as post_creat, (l.charttime - h.he_offset) as last_creat_length
    FROM hypo_cohort_fluids_vp_subjects h
    LEFT JOIN labevents l ON
    H.SUBJECT_ID = L.SUBJECT_ID 
    WHERE l.itemid = 50912 AND  ( l.charttime - h.he_offset) between '00:00:00' and '72:00:00'
    order by SUBJECT_ID
)
, lastcreat AS(
  SELECT h.*
    FROM post72creat h 
    INNER JOIN
    (
        SELECT subject_id, MAX(POST_CREAT) MAXcreat
        FROM post72creat
        GROUP BY subject_id
    ) t ON h.subject_id = t.subject_id AND h.post_creat = t.MAXcreat
)
--select * from lastcreat order by subject_id;
, post72creat_subj as(
    select distinct on (subject_id) * from lastcreat order by subject_id
)
select * 
into hypo_post_creat
from post72creat_subj;
--5,421 patients.

/* Join Two creatinine tables together. */
WITH hypo_trt_creat AS(
    select b.*, p.post_creat_time, p.post_creat, p.last_creat_length
    from hypo_cohort_basecreat b
    LEFT JOIN hypo_post_creat p
    ON b.icustay_id = p.icustay_id
)
select *, (post_creat - pre_creat) as rise_creat 
into hypo_trt_creat
from hypo_trt_creat;

--5,162 patients with complete creatinine data.

/* Create table with all patients */
WITH no_creat AS(
	select * from hypo_cohort_fluids_vp_subjects where subject_id not in (select subject_id from hypo_trt_creat)
)
, all_subj AS(
    select *, NULL AS pre_creat_time, NULL AS pre_creat, NULL AS first_creat_length, NULL AS post_creat_time, NULL AS post_creat, NULL AS last_creat_length, NULL AS rise_creat 
    from no_creat
    union all
    select * from hypo_trt_creat
)
select distinct on (subject_id) * into hypo_trt_creat_all from all_subj;
-----------------------------------------------------------------------------------------------------------------------------------------------

/* Step 2 - SAPS I Info, Service Type */

select * from d_items where lower(label) like '%service%'; -- 1125, 224640

With all_saps as(
select * from sapsi_table where icustay_id in (select icustay_id from hypo_trt_creat_all)
and    age_score is not null and
 	     hr_score is not null and
         sysbp_score is not null and
         resp_score is not null and
         temp_score is not null and
         uo_score is not null and
         vent_score is not null and
         bun_score is not null and
         hematocrit_score is not null and
         wbc_score is not null and
         glucose_score is not null and
         potassium_score is not null and
         sodium_score is not null and
         bicarbonate_score is not null and
         gcs_score is not null
    )
, saps_info as(
	select h.*, saps
    from hypo_trt_creat_all h
    inner join all_saps s
    on h.icustay_id = s.icustay_id
)
select * into hypo_saps from saps_info order by subject_id; -- 5,329

/* Info on SAPS
select k, percentile_disc(k) within group (order by saps)
from saps_info, generate_series(0.25, 0.75, 0.25) as k 
--where dbsource='metavision'
group by k;

select * from sapsi_table;

*/
WITH recservice as (
select distinct h.*, c.value as service,
       coalesce(last_value(c.charttime) over (partition by c.subject_id order by c.charttime 
      rows between unbounded preceding and unbounded following),h.he_onset) reccharttime
        from hypo_saps h
        LEFT join chartevents_adult c
        on h.icustay_id = c.icustay_id
        and c.itemid in (1125,224640)
        and c.charttime < h.he_onset
)
--select * from recservice order by subject_id;

, allservice as(
  select distinct h.subject_id, h.icustay_id, h.he_onset, h.he_offset,
  coalesce(c.value, h.service) as servicemod
  from recservice h left join chartevents_adult c on h.subject_id = c.subject_id and c.itemid in (1125,224640)
     and c.charttime between h.reccharttime and h.he_offset
  order by h.subject_id
)
--select * from allservice;

, servicecount as(
select subject_id, icustay_id, count(*) as scount from allservice group by subject_id, icustay_id
)
--select * from servicecount;

, finalservice as(
select s.subject_id, s.icustay_id, r.he_onset, r.he_offset, r.servicemod from servicecount s,
  allservice r where s.subject_id = r.subject_id and scount < 2 order by subject_id
)
select h.*, f.servicemod 
into hypo_creat_serv
from hypo_saps h
left join finalservice f
on h.icustay_id = f.icustay_id
where f.servicemod is not null;

--------------------------------------------------------------------------------------------------------------------------------------------------------------------

/* Step 3. 

mean MAP in the 3 h period immediately prior to the HE onset (measure of haemodynamic status before the HE), 

total volume of urine output in the 3-h period immediately prior to the HE onset >200 ml (surrogate measure of organ perfusion) 

https://github.com/MIT-LCP/mimic-code has code to find itemid's for urine output.
*/

-- Get mean MAP and total urine output (ml) 3 hours prior to HE onset.
with map_ds as (
    select h.icustay_id, h.he_onset , avg(c.valuenum) as mean_MAP
    from hypo_creat_serv h 
    inner join chartevents_adult c on h.icustay_id = c.icustay_id
    where c.itemid in (220052, 52, 456, 220181)
          and ( h.he_onset - c.charttime ) between '00:00:00' and '03:00:00'
    group by h.icustay_id, h.he_onset
)
--select * from map_ds;
, tuv_ds as (
    select  m.icustay_id, m.he_onset, m.mean_MAP,  SUM(
    -- we consider input of GU irrigant as a negative volume
    case when oe.itemid = 227488 then -1*value
    else value end
    ) as tot_urine_vol
    from outputevents oe
    inner join map_ds m on m.icustay_id = oe.icustay_id
    where oe.itemid in (
  -- these are the most frequently occurring urine output observations in CareVue
  40055, -- "Urine Out Foley"
  43175, -- "Urine ."
  40069, -- "Urine Out Void"
  40094, -- "Urine Out Condom Cath"
  40715, -- "Urine Out Suprapubic"
  40473, -- "Urine Out IleoConduit"
  40085, -- "Urine Out Incontinent"
  40057, -- "Urine Out Rt Nephrostomy"
  40056, -- "Urine Out Lt Nephrostomy"
  40405, -- "Urine Out Other"
  40428, -- "Urine Out Straight Cath"
  40086,--	Urine Out Incontinent
  40096, -- "Urine Out Ureteral Stent #1"
  40651, -- "Urine Out Ureteral Stent #2"

  -- these are the most frequently occurring urine output observations in MetaVision
  226559, -- "Foley"
  226560, -- "Void"
  226561, -- "Condom Cath"
  226584, -- "Ileoconduit"
  226563, -- "Suprapubic"
  226564, -- "R Nephrostomy"
  226565, -- "L Nephrostomy"
  226567, --	Straight Cath
  226557, -- R Ureteral Stent
  226558, -- L Ureteral Stent
  227488, -- GU Irrigant Volume In
  227489  -- GU Irrigant/Urine Volume Out
)
		and oe.value between 0 and 5000 -- sanity check on urine value
        and oe.icustay_id is not null
        and ( m.he_onset - oe.charttime ) between '00:00:00' and '03:00:00'
    group by m.icustay_id, m.he_onset, m.mean_MAP
        )
select h.*, tuv_ds.mean_map, tuv_ds.tot_urine_vol
into he_covariates_ds1
from hypo_creat_serv h
inner join tuv_ds on h.icustay_id = tuv_ds.icustay_id
order by h.subject_id; -- 3,513 subjects

------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/* Step 4  -- Add elixhauser score. */

With hadm_icu_elix as (
	select ELIXHAUSER_QUAN_SCORE.hadm_id, icustays.icustay_id, elixhauser_vanwalraven, elixhauser_sid29, elixhauser_sid30   
    from ELIXHAUSER_QUAN_SCORE
    inner join icustays 
    on ELIXHAUSER_QUAN_SCORE.hadm_id = icustays.hadm_id
    )
select h.*, elixhauser_vanwalraven, elixhauser_sid29, elixhauser_sid30
into HYPO_DS_ELIX_3hr_covs
from he_covariates_ds1 h
inner join hadm_icu_elix k
on k.icustay_id = h.icustay_id;

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/* Step 5 -- Include weight */

select h.*, w.weight
into HYPO_DS_ELIX_3hr_covs2
from HYPO_DS_ELIX_3hr_covs h
LEFT join weight_table w
on w.icustay_id = h.icustay_id
; -- 3,306 subjects have weight information.

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/* Step 6 -- Inlcude heart rate, SPO2. */


-- add more covariates

with hr_ds as (
    select h.icustay_id, h.he_onset, avg(ce.valuenum) as mean_heart_rate
    from HYPO_DS_ELIX_3hr_covs2 h
    inner join chartevents_adult ce on h.icustay_id = ce.icustay_id
    where ce.itemid in (211, 220045 )
          and ( h.he_onset - ce.charttime ) between '00:00:00' and '03:00:00'
    group by h.icustay_id, h.he_onset
)
select h.*, hr_ds.mean_heart_rate
into HYPO_DS_ELIX_3hr_covs3
from HYPO_DS_ELIX_3hr_covs2 h
inner join hr_ds on h.icustay_id = hr_ds.icustay_id
;


with spO2_ds as (
    select h.icustay_id, h.he_onset, avg(ce.valuenum) as mean_spO2
    from HYPO_DS_ELIX_3hr_covs3 h
    inner join chartevents_adult ce on h.icustay_id = ce.icustay_id
    where ce.itemid in ( 646, 220277)
          and ( h.he_onset - ce.charttime ) between '00:00:00' and '03:00:00'
    group by h.icustay_id, h.he_onset
)
select h.*, spO2_ds.mean_spO2
into HYPO_DS_ELIX_3hr_covs4
from HYPO_DS_ELIX_3hr_covs3 h
LEFT join spO2_ds on h.icustay_id = spO2_ds.icustay_id
; -- 3,498 subjects with SPO2 information.

with resp_ds as (
    select h.icustay_id, h.he_onset, avg(ce.valuenum) as mean_resp
    from HYPO_DS_ELIX_3hr_covs3 h
    inner join chartevents_adult ce on h.icustay_id = ce.icustay_id
    where ce.itemid in ( 
    					615,618,220210,224690
    					)
          and ( h.he_onset - ce.charttime ) between '00:00:00' and '03:00:00'
    group by h.icustay_id, h.he_onset
)
select h.*, resp_ds.mean_resp
into HYPO_DS_ELIX_3hr_covs5
from HYPO_DS_ELIX_3hr_covs4 h
LEFT join resp_ds on h.icustay_id = resp_ds.icustay_id
; -- 3,492 subjects with SPO2 information.


--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/* Step 7. Convert age into numeric form. */

select h.*, extract( EPOCH from age) / 60 / 60 / 24 / 365.25 as age_numeric
into HYPO_DS_ELIX_3hr_covs6
from HYPO_DS_ELIX_3hr_covs5 h;

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/* Step 8. Lastly, obtain current service in service table. Provides more detail about the services received by suubject during Hospital Stay. */

WITH curr_service AS(
	select h.*, s.transfertime, s.curr_service,
		coalesce(last_value(s.transfertime) over (partition by s.hadm_id order by s.transfertime 
      	rows between unbounded preceding and unbounded following), s.transfertime) rectrantime
    from HYPO_DS_ELIX_3hr_covs5 h
    inner join services s
    on h.hadm_id = s.hadm_id
    where s.transfertime < h.he_onset
    order by h.subject_id
)
select h.*, s.curr_service
into HYPO_DS_ELIX_3hr_covs7
from HYPO_DS_ELIX_3hr_covs6 h
inner join curr_service s 
on s.hadm_id = h.hadm_id
where transfertime = rectrantime;
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/* Remove subjects with missing baseline creatinine levels. */

select * into HYPO_DS_ELIX_3hr_covs8 from HYPO_DS_ELIX_3hr_covs7 where pre_creat is not null and elixhauser_vanwalraven is not null;

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/* Final Cohort */

select * into HYPO_FINAL_DS from HYPO_DS_ELIX_3hr_covs8 order by subject_id;

--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/* The code below provides simple demographic and clinical characteristics of subjects in our dataset. */

/* Number of subjects from each database. */
select count(*) from HYPO_FINAL_DS
where dbsource != 'metavision'; -- 2,120 subjects from carevue.


select count(*) from HYPO_FINAL_DS
where dbsource = 'metavision'; -- 1,252 subjects from metavision.

select count(*) from HYPO_FINAL_DS; -- 3,372 subjects in entire database.

/* Breakdown by treatment */
-- carevue
select count(*) from HYPO_FINAL_DS
where dbsource != 'metavision' and iv_fluid_ind = 1 and vp_ind = 0 and rise_creat is not null; -- 210 subjects from carevue with IV Fluid administration.


select count(*) from HYPO_FINAL_DS
where dbsource != 'metavision' and iv_fluid_ind = 0 and vp_ind = 1 and rise_creat is not null;  -- 410 subjects from carevue with vasopressor treatment.

-- metavision
select count(*) from HYPO_FINAL_DS
where dbsource = 'metavision' and iv_fluid_ind = 1 and vp_ind = 0 and rise_creat is not null;  -- 252 subjects from metavision with IV Fluid administration.


select count(*) from HYPO_FINAL_DS
where dbsource = 'metavision' and iv_fluid_ind = 0 and vp_ind = 1 and rise_creat is not null;  -- 73 subjects from metavision with vasopressor treatment.

/* Full Data set */
select count(*) from HYPO_FINAL_DS
where iv_fluid_ind=1 and vp_ind=0 and rise_creat is not null; -- 462 with IV fluid administration.

select count(*) from HYPO_FINAL_DS
where iv_fluid_ind=0 and vp_ind=1 and rise_creat is not null; -- 483 with vasoactive therapy.

select count(*) from HYPO_FINAL_DS
where iv_fluid_ind=1 and vp_ind=1 and rise_creat is not null; -- 309 with both vasoactive therapy and IV fluid administration.
--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/* Examine Full data set */

select count(*) from HYPO_FINAL_DS where gender='M'; --1,876
select count(*) from HYPO_FINAL_DS where HOSPITAL_EXPIRE_FLAG=1; --358
select count(*) from HYPO_FINAL_DS; --3372

/* Find quartiles for he_length, LOS, AGE */
select k, percentile_disc(k) within group (order by he_length)
from HYPO_FINAL_DS, generate_series(0.25, 0.75, 0.25) as k 
group by k;

select k, percentile_disc(k) within group (order by los)
from HYPO_FINAL_DS, generate_series(0.25, 0.75, 0.25) as k 
group by k;

select k, percentile_disc(k) within group (order by age_numeric)
from HYPO_FINAL_DS, generate_series(0.25, 0.75, 0.25) as k 
group by k;

select k, percentile_disc(k) within group (order by saps)
from HYPO_FINAL_DS, generate_series(0.25, 0.75, 0.25) as k 
group by k;


/* Service breakdown in Full Data */
select servicemod, count(subject_id) from HYPO_FINAL_DS group by servicemod;

/* Service breakdown by Services table */
select curr_service, count(subject_id) from HYPO_FINAL_DS group by curr_service;


---------------------------------------------------------------------------------------------------------------------------------------------------------

/* Examine Carevue data set */

select count(*) from HYPO_FINAL_DS where dbsource !='metavision' and gender='M'; -- 1,194
select count(*) from HYPO_FINAL_DS where dbsource !='metavision' and HOSPITAL_EXPIRE_FLAG=1; -- 226
select count(*) from HYPO_FINAL_DS where dbsource !='metavision'; --2128

/* Find quartiles for he_length, LOS, AGE */
select k, percentile_disc(k) within group (order by he_length)
from HYPO_FINAL_DS, generate_series(0.25, 0.75, 0.25) as k 
where dbsource !='metavision'
group by k;

select k, percentile_disc(k) within group (order by los)
from HYPO_FINAL_DS, generate_series(0.25, 0.75, 0.25) as k
where dbsource !='metavision'
group by k;

select k, percentile_disc(k) within group (order by age_numeric)
from HYPO_FINAL_DS, generate_series(0.25, 0.75, 0.25) as k 
where dbsource !='metavision'
group by k;

select k, percentile_disc(k) within group (order by saps)
from HYPO_FINAL_DS, generate_series(0.25, 0.75, 0.25) as k 
where dbsource !='metavision'
group by k;


/* Service breakdown in CareVue Data */
select servicemod, count(subject_id) from HYPO_FINAL_DS where dbsource!='metavision' group by servicemod;

/* Service breakdown in CareVue Data by Services table */
select curr_service, count(subject_id) from HYPO_FINAL_DS where dbsource!='metavision' group by curr_service;

---------------------------------------------------------------------------------------------------------------------------------------------------------

/* Examine Metavision data set */

select count(*) from HYPO_FINAL_DS where dbsource ='metavision' and gender='M'; -- 687
select count(*) from HYPO_FINAL_DS where dbsource ='metavision' and HOSPITAL_EXPIRE_FLAG=1; -- 133
select count(*) from HYPO_FINAL_DS where dbsource ='metavision'; --1253

/* Find quartiles for he_length, LOS, AGE */
select k, percentile_disc(k) within group (order by he_length)
from HYPO_FINAL_DS, generate_series(0.25, 0.75, 0.25) as k 
where dbsource ='metavision'
group by k;

select k, percentile_disc(k) within group (order by los)
from HYPO_FINAL_DS, generate_series(0.25, 0.75, 0.25) as k
where dbsource ='metavision'
group by k;

select k, percentile_disc(k) within group (order by age_numeric)
from HYPO_FINAL_DS, generate_series(0.25, 0.75, 0.25) as k 
where dbsource ='metavision'
group by k;

select k, percentile_disc(k) within group (order by saps)
from HYPO_FINAL_DS, generate_series(0.25, 0.75, 0.25) as k 
where dbsource ='metavision'
group by k;

/* Service breakdown in Data */
select servicemod, count(subject_id) from HYPO_FINAL_DS where dbsource='metavision' group by servicemod;

/* Service breakdown in Metavision Data by Services table */
select curr_service, count(subject_id) from HYPO_FINAL_DS where dbsource='metavision' group by curr_service;

---------------------------------------------------------------------------------------------------------------------------------------------------------

/* SAPS I score in carevue*/
select k, percentile_disc(k) within group (order by saps)
from sapsi_table, generate_series(0.25, 0.75, 0.25) as k 
where icustay_id in (select icustay_id from icustays where dbsource!='metavision')
group by k;


/* SAPS I score in metavision */
select k, percentile_disc(k) within group (order by saps)
from sapsi_table, generate_series(0.25, 0.75, 0.25) as k 
where icustay_id in (select icustay_id from icustays where dbsource='metavision')
group by k;

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/* Final table for HE subjects */
select * from HYPO_FINAL_DS;

-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
/* Drop Extra Tables */
drop table HYPO_DS_ELIX_3hr_covs7;
drop table HYPO_DS_ELIX_3hr_covs6;
drop table HYPO_DS_ELIX_3hr_covs5;
drop table HYPO_DS_ELIX_3hr_covs4;
drop table HYPO_DS_ELIX_3hr_covs3;
drop table HYPO_DS_ELIX_3hr_covs2;
drop table HYPO_DS_ELIX_3hr_covs;

---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

/* Full Data set  with available weight. */

select count(*) from HYPO_FINAL_DS
where iv_fluid_ind=1 and vp_ind=0 and rise_creat is not null and weight is not null; -- 436 with IV fluid administration.

select count(*) from HYPO_FINAL_DS
where iv_fluid_ind=0 and vp_ind=1 and rise_creat is not null and weight is not null; -- 467 with vasoactive therapy.

