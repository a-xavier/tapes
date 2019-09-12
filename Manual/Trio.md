# ```--trio```

The ```--trio``` option allows TAPES to assign the PS2 criteria, detecting de-novo mutations from offspring in trios with healthy parents.

A few clarifications on ```--trio``` :  
- ```--trio``` will remove from the final report the "healthy" parents in trios
- You **can** use the same sample in different trios (for example in a case of multiple siblings) but you probably should not use an individual as case offspring and also control parent, as it would mean the individual in not a "healthy" parent.
- In the main report, a new column for each trio will be created using the family ID. **The PS2 criteria will be assigned globally if one of the variant is de-novo in any trio in the cohort**.