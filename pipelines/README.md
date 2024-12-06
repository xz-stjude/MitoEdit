# Editing Pipelines

### 1. Mok2020_G1397

#### Purpose
This pipeline predicts target windows for **C.G-to-T.A** conversions in a **5'-TC context**.

#### Steps
1. Find **TC** and **GA** contexts for the target **C** (top strand) or **G** (bottom strand) in the DNA sequence.
2. For **TC context**, generate target windows (14-18bp) with the target **C** at positions 4-7 from the **3' end** (on the top strand).
3. For **GA context**, generate target windows (14-18bp) with the target **G** at positions 4-7 from the **5' end** (on the top strand).
4. Extract **60bp adjacent sequence** (30bp on each side of the target base).
5. Identify and mark potential **bystander edits** in the windows within the same targeting window.
6. Mark target bases with **[ ]** and bystander bases with **{ }**.

### 2. Mok2020_G1333




### 3. Mok2022_G1397_DddA11



### 4. Cho_sTALEDs
