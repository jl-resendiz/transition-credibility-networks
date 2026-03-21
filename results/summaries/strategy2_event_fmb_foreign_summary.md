# Strategy 2: Event-level Fama-MacBeth and Foreign Exposure

Events (first-mover): 179
Daily obs: 122838; Monthly obs (+12): 76819

## Fama-MacBeth (Daily [-1,+20])
- N_events: 179
- mean b: -0.9814
- SE: 0.3628
- t: -2.70

## Fama-MacBeth (Monthly [-1,+12])
- N_events: 179
- mean b: -447.9870
- SE: 448.4902
- t: -1.00

## Foreign exposure (Daily)
- w_foreign: -2.2450 (SE_event=0.5254, SE_tw=1.7290)
- w_domestic: -0.2751 (SE_event=0.3540, SE_tw=0.5547)
- same_sector: 0.0639 (SE_event=0.0105, SE_tw=0.0450)
- N=122838, R2=0.0004

## Foreign exposure (Monthly)
- w_foreign: -0.2431 (SE_event=0.5992, SE_tw=1.2275)
- w_domestic: -0.1359 (SE_event=0.2241, SE_tw=0.3482)
- same_sector: 0.0246 (SE_event=0.0129, SE_tw=0.0238)
- N=76819, R2=0.0001

