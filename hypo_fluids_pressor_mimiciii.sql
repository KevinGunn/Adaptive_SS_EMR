set search_path to mimiciii;

/* 

First step is to find itemid for Fluid treatment in MIMIC-III.
We need to find itemid for carevue and metavision datbases seperately.

*/


/*

Obtain Treatment patient received between IV fluids and Vasopressors. Also find hypotensive episode duration post treatment.
We include treatment events up to 59 minutes after HE to allow for inclusion of fluid items that were administered in HE window but are entered
after fluid administration is completed, but when charttime may not be within HE window.

*/


/* Find all patients that recieved IV fluid treatments. */

--patients in cv.
WITH fluid_patients_cv AS (
     select h.icustay_id , h.he_onset, h.he_offset, inputevents_cv.itemid as cvfluid_id, inputevents_cv.charttime as cvfluid_input_time, inputevents_cv.amount as cvfluid_vol
     from hypo_cohort_final h
     inner join inputevents_cv on h.icustay_id = inputevents_cv.icustay_id
     where h.dbsource != 'metavision' AND inputevents_cv.charttime BETWEEN h.he_onset AND (h.he_offset + INTERVAL '59 min')
            AND inputevents_cv.itemid in ( 
                  30160 -- D5 Normal Saline 
                , 30014 -- D5/.2NS
                , 30015 -- D5/.45NS
                , 30060 -- D5NS
                , 45374 -- D5NS 40k
                , 45358 -- D%NS with 40 meqkcl
                , 30061 -- D5RL
                ,30176 --0.25% NS
                ,30161 -- 0.3% NS
                , 30020 -- 0.45% NS
                , 30353 -- 0.45% NS
                , 30352 -- 0.9% NS
                , 30143 -- 3% NS
                ,  30018 --	.9% Normal Saline
                  , 30021 --	Lactated Ringers
                  , 40850 --	ns bolus
                  , 41491 --	fluid bolus
                  , 42639 --	bolus
                  , 43819 --	1:1 NS Repletion.
                  , 30063 --	IV Piggyback
                  , 44160 --	BOLUS
                  , 42383 --	cc for cc replace
                  , 42297 --	Fluid bolus
                  , 42453 --	Fluid Bolus
                  , 41490 --	NS bolus
                  , 41371 --	ns fluid bolus
                  , 41322 --	rl bolus
                  , 41896 --	ivf boluses
                  , 41428 --	ns .9% bolus
                  , 44200 --	FLUID BOLUS
                  , 42844 --	NS fluid bolus
                  , 41356 --	IV Bolus
                  , 42548 --	NS Bolus
                  , 44184 --	LR Bolus
                  , 44521 --	LR bolus
                  , 44741 --	NS FLUID BOLUS
                  , 44126 --	fl bolus
                  , 44110 --	RL BOLUS
                  , 44633 --	ns boluses
                  , 44983 --	Bolus NS
                  , 44815 --	LR BOLUS
                  , 43986 --	iv bolus
                  , 45079 --	500 cc ns bolus
                  , 46781 --	lr bolus
                  , 45155 --	ns cc/cc replacement
                  , 41467 --	NS IV bolus
                  , 44367 --	LR
                  , 40423 --	Bolus
                  , 44263 --	fluid bolus ns
                  , 42749 --	fluid bolus NS
                  , 45480 --	500cc ns bolus
                  , 44491 --	.9NS bolus
                  , 41695 --	NS fluid boluses
                  , 41392 --	ns b
                  , 45989 --	NS Fluid Bolus
                  , 45137 --	NS cc/cc
                  , 44053 --	normal saline bolus
                  , 41237 --	ns fluid boluses
                  , 44426 --	bolus ns
                  , 44894 --	N/s 500 ml bolus
                  , 41380 --	nsbolus
 		)
    			AND amount is not NULL
        		AND coalesce(amount,0) > 250
                AND coalesce(amount,0) < 2000
    
    		OR h.dbsource != 'metavision' and inputevents_cv.charttime BETWEEN h.he_onset AND (h.he_offset + INTERVAL '59 min')
    			AND inputevents_cv.itemid in (
                    
                   30008 --	Albumin 5%
                  ,30009 --	Albumin 25%
                  ,42832 --	albumin 12.5%
                  ,40548 --	ALBUMIN
                  ,45403 --	albumin
                  ,44203 --	Albumin 12.5%
                  ,30181 -- Serum Albumin 5%
                  ,46564 -- Albumin
                  ,43237 -- 25% Albumin
                  ,43353 -- Albumin (human) 25%

                  ,30012 --	Hespan
                  ,46313 --	6% Hespan

                  ,42975 --	DEXTRAN DRIP
                  ,42944 --	dextran
                  ,46336 --	10% Dextran 40/D5W
                  ,46729 --	Dextran
                  ,40033 --	DEXTRAN
                  ,45410 --	10% Dextran 40
                  ,30011 -- Dextran 40
                  ,30016 -- Dextrose 10%
                  ,42731 -- Dextran40 10%
                  ,30008 -- Albumin 5%
                  ,30016 -- Dextrose 10%
                  
                    )
   ),
 --select distinct on (icustay_id) * from fluid_patients_cv;
CV_FLUID_ID AS(
	select distinct icustay_id from fluid_patients_cv
)
--select icustay_id into from CV_FLUID_ID;

/* Patients can be given treatment multiple times so by using DISTINCT we find each patient that recieved treatment once. */


, fluid_pats_ce AS (
    select h.icustay_id , ce.itemid , ce.charttime as cvfluid_input_time, ce.valuenum as cvfluid_vol
     from hypo_cohort_final h
     inner join chartevents_adult ce on h.icustay_id = ce.icustay_id
     where h.dbsource != 'metavision' and ce.charttime BETWEEN h.he_onset AND (h.he_offset + INTERVAL '59 min') 
    AND ce.itemid in
  ( 
      2510 --	DEXTRAN LML 10%
    , 3087 --	DEXTRAN 40  10%
    , 6937 --	Dextran
    , 3087 -- | DEXTRAN 40  10%
    , 3088 --	DEXTRAN 40%
      -- NOTE - THESE ARE COLLOIDS IN CE TABLE.
  )
  --and ce.valuenum is not null
    OR h.dbsource != 'metavision' and ce.charttime BETWEEN h.he_onset AND (h.he_offset + INTERVAL '59 min') 
    AND ce.itemid in
    (
      
     1460 -- D5W 100meq NaHCO3
    ,5392 -- D10W bolus
    ,1848 -- D5W (Dextrose in water)
    ,2619 --3% NS
    ,5321 -- ns bolus
    ,1634 -- lactated ringers
    ,4647 -- NS bolus
    ,4970 -- NS bolus
    ,5199 --NS bolus
    ,5333 -- saline flush
    ,6190 --NS
    )
    AND valuenum > 250
),

CV_FLUID_ID_CE AS(
	select distinct icustay_id from fluid_pats_ce
)

-- patients in mv.    
, fluid_patients_mv AS (
       select h.icustay_id , inputevents_mv.itemid as mvfluid_id, inputevents_mv.starttime as mvfluid_input_time, inputevents_mv.amount as mvfluid_vol ,
    		   inputevents_mv.cancelreason
            from hypo_cohort_final h
        	inner join inputevents_mv on h.icustay_id = inputevents_mv.icustay_id
        		where h.dbsource = 'metavision' and inputevents_mv.starttime BETWEEN h.he_onset AND (h.he_offset + INTERVAL '59 min')
                	AND inputevents_mv.itemid in ( 
                    	225158, -- NaCl 0.9%
    					225828, -- LR
    					225944, -- Sterile Water
    					225797,  -- Free Water
                        225823, -- D5 1/2NS
                        225825, -- D5NS
                        225827, -- D5LR
                        225941 -- D5 1/4NS
                    )
						
        			AND
                      -- in MetaVision, these ITEMIDs appear with a null rate IFF endtime=starttime + 1 minute
                      -- so it is sufficient to:
                      --    (1) check the rate is > 240 if it exists or
                      --    (2) ensure the rate is null and amount > 240 ml
					    (
                          (rate is not null and rateuom = 'mL/hour' and rate > 250)
                          OR (rate is not null and rateuom = 'mL/min' and rate > (250/60.0))
                          OR (rate is null and amountuom = 'L' and amount > 0.250)
                          OR (rate is null and amountuom = 'ml' and amount > 250)
                        )
    				AND cancelreason = 0
                    AND statusdescription != 'Rewritten'
    			OR h.dbsource = 'metavision' AND inputevents_mv.starttime BETWEEN h.he_onset AND (h.he_offset + INTERVAL '59 min')
                	AND inputevents_mv.itemid in ( 
                    220864, --	Albumin 5%	7466 132 7466
                    220862, --	Albumin 25%	9851 174 9851
                    225174, --	Hetastarch (Hespan) 6%	82 1 82
                    225795,  --	Dextran 40	38 3 38
                    225796 --  Dextran 70
                    -- below ITEMIDs not in use
                   -- 220861 | Albumin (Human) 20%
                   -- 220863 | Albumin (Human) 4%
                   ,220949 -- Dextrose 5%
                   ,220950 -- Dextrose 10%
                   ,220952 -- Dextrose 50%
                    )
                    
    			    AND cancelreason = 0
    				AND statusdescription != 'Rewritten'
    				AND rate is not null
    				AND rate > 0
        ),
icustays_mv_fluid as (
	select distinct icustay_id from fluid_patients_mv
),
fluid_icu_ids as(
    select icustay_id from CV_FLUID_ID
    union
    select * from icustays_mv_fluid
    union 
    select * from CV_FLUID_ID_CE
)
select distinct icustay_id into fluid_icu_ids 
from fluid_icu_ids;

-- now make indicator variable if patient recieved fluids treatment.
select *, ( CASE
		WHEN icustay_id in (select icustay_id from fluid_icu_ids ) Then 1
    	ELSE 0
    	END ) AS IV_fluid_ind
into hypo_cohort_fluids_trt
from hypo_cohort_final
order by subject_id;

------	------	------	------	------	------	------	------	------	-----	-----	-----	-----	-----
                
/* Next, examine patients who recieved vassopressor treatment */ 


--patients in cv.
WITH vp_patients_cv AS (
     select h.icustay_id , h.he_onset, h.he_offset, ii.itemid , ii.charttime as cvvp_input_time, ii.rate as cvvp_rate
     from hypo_cohort_final h
     inner join inputevents_cv ii on h.icustay_id = ii.icustay_id
     where h.dbsource != 'metavision' and ii.charttime BETWEEN h.he_onset AND h.he_offset --+ INTERVAL '59 min')
            AND  ii.itemid in
				(
                   30047,30120 -- norepinephrine
    			  ,30044,30119,30309 -- epinephrine 
    			  ,30127,30128 -- phenylephrine
    		      ,30051, 42273, 42802 -- vasopressin -- ,42273, 46570 -- vasopressin
    		      ,30306, 30042 -- dobutamine
                  ,30043,30307,30125 -- dopamine
					)
    		 AND rate is not null
  			 AND rate > 0
        ),
        --select distinct on (icustay_id) * from vp_patients_cv; 
icustays_cv_vp_ie as (
	select distinct icustay_id from vp_patients_cv 
)

-- patients in mv.    
, vp_patients_mv AS (
        select h.icustay_id , inputevents_mv.itemid, inputevents_mv.starttime , inputevents_mv.amount as mvvp_amount ,
    		   inputevents_mv.cancelreason
            from hypo_cohort_final h
            inner join inputevents_mv on h.icustay_id = inputevents_mv.icustay_id
        	WHERE h.dbsource = 'metavision' and inputevents_mv.starttime BETWEEN h.he_onset AND (h.he_offset + INTERVAL '59 min') 
                AND inputevents_mv.itemid in
				(
                    221906 -- norepinephrine
  					,221289 -- epinephrine
  					,221749 -- phenylephrine
  					,222315 -- vasopressin
  					,221653 -- dobutamine
                    ,221662 -- dopamine
                    ,222315 -- vasopressin
                    --,221986 -- milrinone
                    --,227692 -- isuprel
					)
    			AND rate is not null
  				AND rate > 0
  				AND statusdescription != 'Rewritten' -- only valid orders
    			AND cancelreason = 0
        ), 
icustays_mv_vp as (
	select distinct icustay_id from vp_patients_mv
), 
vp_icu_ids as(
    select * from icustays_mv_vp
    union
    select * from icustays_cv_vp_ie
)
select * into vp_icu_ids from vp_icu_ids;

-- now make indicator variable if patient recieved fluids treatment.
select *, ( CASE
		WHEN icustay_id in (select icustay_id from vp_icu_ids ) Then 1
    	ELSE 0
    	END ) AS vp_ind
into hypo_cohort_fluids_vp
from hypo_cohort_fluids_trt;

select * from hypo_cohort_fluids_vp;

-- Final Table. "Distinct on Subject_id" removes Subject Duplicates i.e. subjects that belong to both Metavision and Carevue.
SELECT distinct on(subject_id) * INTO hypo_cohort_fluids_vp_subjects FROM hypo_cohort_fluids_vp order by subject_id;



-- Breakdown of subjects belonging to each cohort.

-- Subjects who received both in carevue.
select count(*) from hypo_cohort_fluids_vp_subjects where dbsource != 'metavision' and vp_ind=1 and iv_fluid_ind=1; -- 269

-- Subjects who received vasopressors in carevue.
select count(*) from hypo_cohort_fluids_vp_subjects where dbsource != 'metavision' and vp_ind=1 and iv_fluid_ind=0; -- 590

-- Subjects who received IV Fluids in carevue.
select count(*) from hypo_cohort_fluids_vp_subjects where dbsource != 'metavision' and vp_ind=0 and iv_fluid_ind=1; -- 354

-- Subjects who received no treatment in carevue.
select count(*) from hypo_cohort_fluids_vp_subjects where dbsource != 'metavision' and vp_ind=0 and iv_fluid_ind=0; -- 1895

----	------	------

-- Subjects who received both in metavision.
select count(*) from hypo_cohort_fluids_vp_subjects where dbsource = 'metavision' and vp_ind=1 and iv_fluid_ind=1; -- 359

-- Subjects who received vasopressors in metavision.
select count(*) from hypo_cohort_fluids_vp_subjects where dbsource = 'metavision' and vp_ind=1 and iv_fluid_ind=0; -- 164

-- Subjects who received IV Fluids in metavision.
select count(*) from hypo_cohort_fluids_vp_subjects where dbsource = 'metavision' and vp_ind=0 and iv_fluid_ind=1; -- 529

-- Subjects who received no treatment in metavision.
select count(*) from hypo_cohort_fluids_vp_subjects where dbsource = 'metavision' and vp_ind=0 and iv_fluid_ind=0; -- 1545

-----	------	------

-- Subjects who received both.
select count(*) from hypo_cohort_fluids_vp_subjects where vp_ind=1 and iv_fluid_ind=1; -- 628

-- Subjects who received vasopressors.
select count(*) from hypo_cohort_fluids_vp_subjects where vp_ind=1 and iv_fluid_ind=0; -- 754

-- Subjects who received IV Fluids.
select count(*) from hypo_cohort_fluids_vp_subjects where vp_ind=0 and iv_fluid_ind=1; -- 883

-- Subjects who received no treatment
select count(*) from hypo_cohort_fluids_vp_subjects where vp_ind=0 and iv_fluid_ind=0; -- 3440



              
                
                
                
                
                