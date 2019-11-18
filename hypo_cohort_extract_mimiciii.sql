set search_path to mimiciii;

/* 

  Thank you to Dr. Rithi Koshari for sharing some SQL code for extracting patients with hypotensive episodes in MIMIC-II.
  
  This query returns DBSOURCE, HADM_ID, SUBJECT_ID, ICUSTAY_ID, HE_ONSET, HE_OFFSET, HE_LENGTH, LOS, AGE for 
  subjects in MIMIC-III DB on criteria:
  
  Entry: Time when first of 2 measurements of less than 60 mmHg within 1 hour was made
  Exit: Time when first of 2 measurements of greater than 60 mmHg within 1 hour was made   
    
****
  In this query, we find subjects who fit both entry and exit criteria and experienced one hypotensive episode duirng ICU Stay.
  Subjects with one hypotensive episode during their ICU Stay are included in our cohort. 
****
  
Hypoentries: gets charttimes and values that qualify as entry criteria
  
Hypoentries2: This table finds all patients that fit the definition of the start of a hypotensive episode.

Hypoexits: gets charttimes and values that qualify as exit criteria

Hypoexits2: This table finds all patients that fit the definition of the end of a hypotensive episode.


Allevents: Gathers all ENTRY and EXIT events into one table, as well as
  information about the previous record for filtering
  
Hypoevents: Find all patients that expierenced a HE.
   
HEtable1/final query: Assembles an HE event. Contains all ENTRY times and EXIT times.
  
After this query, the SQL file "Hypo_fluid_pressor_mimiciii.sql" will use hypo_cohort_final to obtain patients who recieved treatment.

*/


/* 

This code finds all adult patients and makes chartevents tabe with their information. 
I use this table pull data from chartevents_adult. The code below creates this table.

/* Remove patients younger than 15. */
WITH first_admission_time AS (
SELECT 
    p.subject_id, p.dob, p.gender, 
    MIN (a.admittime) AS first_admittime
FROM patients p
INNER JOIN admissions a
ON p.subject_id = a.subject_id
GROUP BY p.subject_id, p.dob, p.gender, a.hadm_id
ORDER BY a.hadm_id, p.subject_id
),
age AS (
SELECT 
    subject_id, dob, gender, first_admittime, 
    age(first_admittime, dob) 
        AS first_admit_age, 
    CASE
       WHEN age(first_admittime, dob) > '89 years'
            THEN '>89'
        WHEN age(first_admittime, dob) >= '15 years'
            THEN 'adult'
        WHEN age(first_admittime, dob) <= '1  year'
            THEN 'neonate'
        ELSE 'middle'
        END AS age_group
FROM first_admission_time
ORDER BY subject_id
)
	,sub_age AS (
    SELECT *
    FROM age
     )
select sub_age.*, icustays.icustay_id
into subject_age
from sub_age
inner join icustays on sub_age.subject_id = icustays.subject_id;
*/

/*
select *
into chartevents_adult
from chartevents
where chartevents.subject_id in (
  select subject_age.subject_id
  from subject_age
  where age_group = 'adult'
	) 
  and chartevents.icustay_id in (
  select icustay_id
  from icustays
  where first_careunit in ('MICU','CCU','SICU','CSRU') or last_careunit in ('MICU','CCU','SICU','CSRU')  
);

*/



CREATE TABLE dual        
	(
     enter_hypo_threshold int, 
     enter_windownum int,
     exit_hypo_threshold int, 
     exit_windownum int
);
INSERT INTO dual
Values (60, 1, 60, 1);

/* BP measurements lower than 40 can be too low to be realistic sometimes. */

With hypoentries as (
select subject_id, icustay_id, itemid, valuenum, charttime as he_time,
    row_number()  over (partition by subject_id, icustay_id order by charttime) as rn,
    last_value(charttime) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowendct,
    last_value(valuenum) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend,
    last_value(itemid) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend_itemid
    FROM chartevents_adult WHERE itemid in (52, 456)
),
--select * from hypoentries limit 50;

hypoentries2 as (
select h.*,
 cast('ENTRY' AS text) as status
 FROM hypoentries h WHERE valuenum is not null and valuenum between 40 and (select enter_hypo_threshold from dual) 
                                      and windowend between 40 and (select enter_hypo_threshold from dual) and icustay_id is not null
    								  and itemid = windowend_itemid and he_time != windowendct
    								  and (windowendct - he_time) <= '01:00:00'
),
--select * from hypoentries2 limit 500;

/* BP measurements greater than 180 can be too low to be realistic sometimes. */

hypoexits as(
	select subject_id, icustay_id, itemid, valuenum, charttime as he_time,
    row_number()  over (partition by subject_id, icustay_id order by charttime) as rn,
    last_value(charttime) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowendct,
    last_value(valuenum) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend,
    last_value(itemid) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend_itemid
    FROM chartevents_adult WHERE itemid in (52, 456) 
),

hypoexits2 as(
    select h.*,
    cast('EXIT' AS text) as status
    FROM hypoexits h WHERE valuenum is not null and valuenum > (select exit_hypo_threshold from dual) and valuenum <= 180
                                       and windowend > (select exit_hypo_threshold from dual) and windowend <= 180 
    								   and icustay_id is not null
    								   and itemid = windowend_itemid and he_time != windowendct
    								   and (windowendct - he_time) <= '01:00:00'
),
--select * from hypoexits2 limit 500;

allevents as(
  select aa.* from (select * from hypoentries2 union select * from hypoexits2) aa order by subject_id, he_time
)
--select * from allevents limit 10000;

, allevents2 as(
	select *,
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lag(itemid,1 ) over (partition by subject_id, icustay_id order by he_time) as prev_itemid
    from allevents
)
--select * from allevents2 limit 10000;
, allevents3_icu as(
        select icustay_id from allevents2 where ( icustay_id = prev_icustay and status != prev_status and itemid = prev_itemid ) 
    								
)
, hypo_events as(
	select * from allevents where icustay_id in (select * from allevents3_icu)
 )
--select * from hypo_events limit 5000;    

, allevents35 as(
  select subject_id, icustay_id, itemid, valuenum, he_time, status, 
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lead(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as next_icustay,
    lead(status, 1) over (partition by subject_id, icustay_id order by he_time) as next_status,
    LEAD(he_time, 1) over (partition by SUBJECT_ID, ICUSTAY_ID order by he_time) as NEXT_CT
  from hypo_events
)
--select * from allevents35 limit 10000; 

, hetable1 as(
  select SUBJECT_ID, ICUSTAY_ID, itemid,
    case 
      when (status = 'ENTRY') 
        then he_time else null 
    end he_onset,
    case 
      when (status = 'ENTRY' and ICUSTAY_ID = NEXT_ICUstay and NEXT_STATUS = 'EXIT')
        then NEXT_CT else null 
    end he_offset
  from allevents35
)
select h.* into he_set 
from hetable1 h where he_onset is not null order by h.subject_id, h.he_onset;

/*
The next block of code finds subjects with one hypotensive episode and retrieves the length of the episode.
*/

with usable_he AS(
SELECT distinct h.subject_id, icustay_id, 
       dense_rank() over (partition by h.subject_id order by h.icustay_id)  + dense_rank() over (partition by h.subject_id order by h.icustay_id desc) -1 as count_icu,
       dense_rank() over (partition by h.subject_id, h.icustay_id order by h.he_onset)  +  dense_rank() over (partition by h.subject_id, h.icustay_id order by h.he_onset desc) -1 
       as count_events,
       dense_rank() over (partition by h.subject_id, h.icustay_id order by h.itemid)  +  dense_rank() over (partition by h.subject_id, h.icustay_id order by h.itemid desc) -1 
       as count_itemid,
       first_value(h.he_offset) OVER (PARTITION BY h.subject_id ORDER BY h.he_offset) first_offset,
       first_value(h.itemid) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) first_item
       FROM he_set h where he_offset is not null
)
--select * from usable_he order by subject_id;
, he_subjs AS(
    select h.* from he_set h, usable_he u where (count_events = 1 and count_icu = 1) 
    and h.he_offset=u.first_offset 
    and he_offset is not null 
    order by subject_id
)
, usable_he2 AS(
SELECT distinct h.subject_id, icustay_id, itemid,
    	 first_value(h.he_onset) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY h.he_onset ) he_onset,
         first_value( h.he_offset ) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY case when h.he_offset is not null then 0 else 1 end ASC, h.he_offset 
         rows unbounded preceding ) he_offset
    FROM he_subjs h where he_offset is not null
)
select *, (he_offset - he_onset) as he_length 
into hypo_cohort_icu 
from usable_he2
order by subject_id;

/*

This block of code removes patients with CMO with 24 hrs of HE. It also attaches dbsource, Hospital Admission ID, ICU intime, LOS, and AGE.

*/


-- Remove pateints with CMO's.
WITH icustay_CMO AS (
    select h.icustay_id, chartevents_adult.itemid, chartevents_adult.value, chartevents_adult.charttime as CMO_time, h.he_length, h.he_onset, h.he_offset 
         from chartevents_adult
         inner join hypo_cohort_icu h on chartevents_adult.icustay_id = h.icustay_id
          where  value in ('Comfort Measures','Comfort measures only','CMO') and
                        (h.he_onset - chartevents_adult.charttime ) between '-24:00:00' and '24:00:00' or
                   value in ('Comfort Measures', 'Comfort measures only','CMO') and
                        (h.he_offset - chartevents_adult.charttime ) between '-24:00:00' and '24:00:00'
    		)
select * into hypo_cohort_final_cv 
from hypo_cohort_icu
where icustay_id not in (select icustay_id from icustay_CMO) order by he_length;

-- intime is admittime to icu for each patient.
with icu_hypo AS(
    select dbsource, hadm_id, h.*, intime, los
        from hypo_cohort_final_cv  h
        INNER JOIN icustays
        ON h.icustay_id = icustays.icustay_id 
        order by h.icustay_id
)
, pats_hypo AS(
 select h.*, gender, dob 
 from icu_hypo h
 INNER JOIN patients p
 ON p.subject_id = h.subject_id
 order by h.subject_id
)
, inhm_hypo AS(
  select p.*, HOSPITAL_EXPIRE_FLAG
  from pats_hypo p
  inner join admissions a
  ON p.hadm_id = a.hadm_id
)
select *, age(he_offset , dob) as age 
into hypo_cohort_cv 
from inhm_hypo;


------------------------------------------------------------------------------------------------------------------------------------------------------
/* The following code repeats the steps above for the patients in the metavision DB. They have different itemids. */
------------------------------------------------------------------------------------------------------------------------------------------------------

/* BP measurements lower than 40 can be too low to be realistic sometimes. */

With hypoentries as (
select subject_id, icustay_id, itemid, valuenum, charttime as he_time,
    row_number()  over (partition by subject_id, icustay_id order by charttime) as rn,
    last_value(charttime) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowendct,
    last_value(valuenum) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend,
    last_value(itemid) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend_itemid
    FROM chartevents_adult WHERE itemid in (220052, 220181)
),
--select * from hypoentries limit 50;

hypoentries2 as (
select h.*,
 cast('ENTRY' AS text) as status
 FROM hypoentries h WHERE valuenum is not null and valuenum between 40 and (select enter_hypo_threshold from dual) 
                                      and windowend between 40 and (select enter_hypo_threshold from dual) and icustay_id is not null
    								  and itemid = windowend_itemid and he_time != windowendct
    								  and (windowendct - he_time) <= '01:00:00'
),
--select * from hypoentries2 limit 500;

/* BP measurements greater than 180 can be too low to be realistic sometimes. */

hypoexits as(
	select subject_id, icustay_id, itemid, valuenum, charttime as he_time,
    row_number()  over (partition by subject_id, icustay_id order by charttime) as rn,
    last_value(charttime) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowendct,
    last_value(valuenum) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend,
    last_value(itemid) over (partition by subject_id, icustay_id order by charttime ROWS BETWEEN CURRENT ROW and 1 following ) as windowend_itemid
    FROM chartevents_adult WHERE itemid in (220052, 220181) 
),

hypoexits2 as(
    select h.*,
    cast('EXIT' AS text) as status
    FROM hypoexits h WHERE valuenum is not null and valuenum > (select exit_hypo_threshold from dual) and valuenum <= 180
                                       and windowend > (select exit_hypo_threshold from dual) and windowend <= 180 
    								   and icustay_id is not null
    								   and itemid = windowend_itemid and he_time != windowendct
    								   and (windowendct - he_time) <= '01:00:00'
),
--select * from hypoexits2 limit 500;

allevents as(
  select aa.* from (select * from hypoentries2 union select * from hypoexits2) aa order by subject_id, he_time
)
--select * from allevents limit 10000;

, allevents2 as(
	select *,
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lag(itemid,1 ) over (partition by subject_id, icustay_id order by he_time) as prev_itemid
    from allevents
)
--select * from allevents2 limit 10000;
, allevents3_icu as(
        select icustay_id from allevents2 where ( icustay_id = prev_icustay and status != prev_status and itemid = prev_itemid ) 
    								
)
, hypo_events as(
	select * from allevents where icustay_id in (select * from allevents3_icu)
 )
--select * from hypo_events limit 5000;    

, allevents35 as(
  select subject_id, icustay_id, itemid, valuenum, he_time, status, 
    lag(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as prev_icustay,
    lag(status, 1) over (partition by subject_id, icustay_id order by he_time) as prev_status,
    lead(icustay_id, 1) over (partition by subject_id, icustay_id order by he_time) as next_icustay,
    lead(status, 1) over (partition by subject_id, icustay_id order by he_time) as next_status,
    LEAD(he_time, 1) over (partition by SUBJECT_ID, ICUSTAY_ID order by he_time) as NEXT_CT
  from hypo_events
)
--select * from allevents35 limit 10000; 

, hetable1 as(
  select SUBJECT_ID, ICUSTAY_ID, itemid,
    case 
      when (status = 'ENTRY') 
        then he_time else null 
    end he_onset,
    case 
      when (status = 'ENTRY' and ICUSTAY_ID = NEXT_ICUstay and NEXT_STATUS = 'EXIT')
        then NEXT_CT else null 
    end he_offset
  from allevents35
)
select h.* 
into he_set2
from hetable1 h where he_onset is not null order by h.subject_id, h.he_onset;

/*
The next block of code finds subjects with one hypotensive episode and retrieves the length of the episode.
*/

with usable_he AS(
SELECT distinct h.subject_id, icustay_id, 
       dense_rank() over (partition by h.subject_id order by h.icustay_id)  + dense_rank() over (partition by h.subject_id order by h.icustay_id desc) -1 as count_icu,
       dense_rank() over (partition by h.subject_id, h.icustay_id order by h.he_onset)  +  dense_rank() over (partition by h.subject_id, h.icustay_id order by h.he_onset desc) -1 
       as count_events,
       dense_rank() over (partition by h.subject_id, h.icustay_id order by h.itemid)  +  dense_rank() over (partition by h.subject_id, h.icustay_id order by h.itemid desc) -1 
       as count_itemid,
       first_value(h.he_offset) OVER (PARTITION BY h.subject_id ORDER BY h.he_offset) first_offset,
       first_value(h.itemid) OVER (PARTITION BY h.subject_id ORDER BY h.he_onset) first_item
       FROM he_set2 h where he_offset is not null
)
--select * from usable_he order by subject_id;
, he_subjs AS(
    select h.* from he_set2 h, usable_he u where (count_events = 1 and count_icu = 1) 
    and h.he_offset=u.first_offset 
    and he_offset is not null 
    order by subject_id
)
, usable_he2 AS(
SELECT distinct h.subject_id, icustay_id, itemid,
    	 first_value(h.he_onset) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY h.he_onset ) he_onset,
         first_value( h.he_offset ) OVER (PARTITION BY h.subject_id, h.icustay_id ORDER BY case when h.he_offset is not null then 0 else 1 end ASC, h.he_offset 
         rows unbounded preceding ) he_offset
    FROM he_subjs h where he_offset is not null
)
select *, (he_offset - he_onset) as he_length 
into hypo_cohort_icu2 
from usable_he2
order by subject_id;

/*

This block of code removes patients with CMO with 24 hrs of HE. It also attaches dbsource, Hospital Admission ID, ICU intime, LOS, and AGE.

*/


-- Remove pateints with CMO's.
WITH icustay_CMO AS (
    select h.icustay_id, chartevents_adult.itemid, chartevents_adult.value, chartevents_adult.charttime as CMO_time, h.he_length, h.he_onset, h.he_offset 
         from chartevents_adult
         inner join hypo_cohort_icu h on chartevents_adult.icustay_id = h.icustay_id
          where  value in ('Comfort Measures','Comfort measures only','CMO') and
                        (h.he_onset - chartevents_adult.charttime ) between '-24:00:00' and '24:00:00' or
                   value in ('Comfort Measures', 'Comfort measures only','CMO') and
                        (h.he_offset - chartevents_adult.charttime ) between '-24:00:00' and '24:00:00'
    		)
select * into hypo_cohort_final_mv 
from hypo_cohort_icu2
where icustay_id not in (select icustay_id from icustay_CMO) order by he_length;

-- intime is admittime to icu for each patient.
with icu_hypo AS(
    select dbsource, hadm_id, h.*, intime, los
        from hypo_cohort_final_mv  h
        INNER JOIN icustays
        ON h.icustay_id = icustays.icustay_id 
        order by h.icustay_id
)
, pats_hypo AS(
 select h.*, gender, dob 
 from icu_hypo h
 INNER JOIN patients p
 ON p.subject_id = h.subject_id
 order by h.subject_id
)
, inhm_hypo AS(
  select p.*, HOSPITAL_EXPIRE_FLAG
  from pats_hypo p
  inner join admissions a
  ON p.hadm_id = a.hadm_id
)
select *, age(he_offset , dob) as age 
into hypo_cohort_mv 
from inhm_hypo;


/* Create final dataset by merging these two datasets. */
with hypo_final AS(
    select * from hypo_cohort_cv
    union
    select * from hypo_cohort_mv
)
select *
into hypo_cohort_final
from hypo_final;


/* Examine data set */

select count(*) from hypo_cohort_cv where gender='M'; --1744
select count(*) from hypo_cohort_cv where HOSPITAL_EXPIRE_FLAG=1; --384
select count(*) from hypo_cohort_cv; --3122

select count(*) from hypo_cohort_mv where gender='M'; --1423
select count(*) from hypo_cohort_mv where HOSPITAL_EXPIRE_FLAG=1; --315
select count(*) from hypo_cohort_mv; --2629

/* Find quartiles for he_length, LOS, AGE */
select k, percentile_disc(k) within group (order by he_length)
from hypo_cohort_final, generate_series(0.25, 0.75, 0.25) as k 
group by k;

select k, percentile_disc(k) within group (order by los)
from hypo_cohort_final, generate_series(0.25, 0.75, 0.25) as k 
group by k;

select k, percentile_disc(k) within group (order by age)
from hypo_cohort_final, generate_series(0.25, 0.75, 0.25) as k 
group by k;



