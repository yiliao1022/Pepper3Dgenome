## Step1: genome alignment with minimap2

#!/bin/bash
~/Softwares/minimap2/minimap2 -c --cs=long ../../ref/CA59 ../../HQ > CA59.HQ.long.paf
~/Softwares/minimap2/minimap2 -c --cs=long ../../ref/CA59 ../../RH89A > CA59.RH89A.long.paf
~/Softwares/minimap2/minimap2 -c --cs=long ../../ref/CA59 ../../RH89B > CA59.RH89B.long.paf
~/Softwares/minimap2/minimap2 -c --cs=long ../../ref/CA59 ../../SL4 > CA59.SL4.long.paf 
~/Softwares/minimap2/minimap2 -c --cs=long ../../RH89A ../../ref/CA59 > RH89A.CA59.long.paf
~/Softwares/minimap2/minimap2 -c --cs=long ../../RH89A ../../HQ > RH89A.HQ.long.paf
~/Softwares/minimap2/minimap2 -c --cs=long ../../RH89A ../../SL4 > RH89A.SL4.long.paf
~/Softwares/minimap2/minimap2 -c --cs=long ../../SL4 ../../ref/CA59 > SL4.CA59.long.paf
~/Softwares/minimap2/minimap2 -c --cs=long ../../SL4 ../../HQ > SL4.HQ.long.paf
~/Softwares/minimap2/minimap2 -c --cs=long ../../SL4 ../../RH89A > SL4.RH89A.long.paf
~/Softwares/minimap2/minimap2 -c --cs=long ../../SL4 ../../RH89B > SL4.RH89B.long.paf

## Step2:
