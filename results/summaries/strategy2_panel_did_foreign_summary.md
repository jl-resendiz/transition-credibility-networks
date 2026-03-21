# Strategy 2 Panel DiD (Foreign Exposure Only)

Event scope: first_mover
Window: tau=[-6,+12], post=[0,+12], overlap=nearest
Exposure: foreign-only (same-country weights set to 0)

Results:
- Event-clustered: beta=+0.0607, SE=0.0615, t=0.99
- Firm-clustered:  beta=+0.0607, SE=0.0217, t=2.79
- Two-way:         beta=+0.0607, SE=0.0135, t=4.48
- N=18,515

Note: two-way SEs are smaller than event-clustered in this stacked design; interpret with caution.
