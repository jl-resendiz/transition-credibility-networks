"""Add announcement_date column to coal_retirement_events.csv.

Markets react to announcements, not physical closures. For a finance event
study, the proper event date is when the retirement was first publicly
disclosed (utility filing, press release, regulatory settlement, board vote).

Announcement dates are loaded from data/raw/events/manual_announcement_dates.csv,
a researcher-curated reference file with columns:
  plant_name, date, source

Sources include:
  - EIA Form 860 planned retirement filings
  - PJM deactivation notices
  - EPA Clean Air Act consent decrees and settlements
  - Utility Integrated Resource Plans and PSC filings
  - Company press releases and board resolutions
  - National energy policy announcements (Energieakkoord, CACJ Act, etc.)

Each date is tagged with source for reproducibility.
"""
import csv, os
import datetime

from _paths import raw_path, derived_path

# ── Load announcement date lookup from CSV ────────────────────────────
# date = YYYY-MM-DD when market first learned of the retirement plan
# For multi-unit plants, use earliest public disclosure.
# Day-15 = month known but not exact day.

ANNOUNCEMENT_DATES = {}
ann_csv = raw_path('events', 'manual_announcement_dates.csv')
with open(ann_csv, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        ANNOUNCEMENT_DATES[row['plant_name']] = (
            row['date'], row['source']
        )


# ── Update the CSV ──────────────────────────────────────────────────
events_path = derived_path('events', 'coal_retirement_events.csv')

all_events = []
with open(events_path, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    fieldnames = list(reader.fieldnames)
    for row in reader:
        all_events.append(row)

# Add announcement_date column if not present
if 'announcement_date' not in fieldnames:
    fieldnames.append('announcement_date')
if 'announcement_source' not in fieldnames:
    fieldnames.append('announcement_source')

updated_count = 0
already_had = 0
not_found = 0

for row in all_events:
    if row['is_first_mover'] != 'True' or row['is_matched'] != 'True':
        row.setdefault('announcement_date', '')
        row.setdefault('announcement_source', '')
        continue

    plant = row['plant_name']

    # Already has an announcement date
    if row.get('announcement_date') and row['announcement_date'].strip():
        already_had += 1
        continue

    if plant in ANNOUNCEMENT_DATES:
        date, source = ANNOUNCEMENT_DATES[plant]
        row['announcement_date'] = date
        row['announcement_source'] = source
        updated_count += 1
    else:
        row.setdefault('announcement_date', '')
        row.setdefault('announcement_source', '')
        not_found += 1

# Write back
with open(events_path, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(all_events)

# Stats
total_fm = sum(1 for r in all_events
               if r['is_first_mover'] == 'True' and r['is_matched'] == 'True')
with_ann = sum(1 for r in all_events
               if r['is_first_mover'] == 'True' and r['is_matched'] == 'True'
               and r.get('announcement_date') and r['announcement_date'].strip())

# Precision: count firm-sourced vs approximate dates
firm_dates = 0
approx_dates = 0
for plant, (d, s) in ANNOUNCEMENT_DATES.items():
    if d.endswith('-15') and 'press release' not in s.lower() and 'settlement' not in s.lower():
        approx_dates += 1
    else:
        firm_dates += 1

# Lead time analysis
import datetime
leads = []
for row in all_events:
    if (row['is_first_mover'] == 'True' and row['is_matched'] == 'True'
            and row.get('announcement_date') and row['announcement_date'].strip()
            and row.get('event_date') and row['event_date'].strip()):
        try:
            ann = datetime.date.fromisoformat(row['announcement_date'])
            ret = datetime.date.fromisoformat(row['event_date'])
            lead_days = (ret - ann).days
            leads.append((row['plant_name'], lead_days))
        except ValueError:
            pass

print(f'=== ANNOUNCEMENT DATE UPDATE ===')
print(f'Total matched first-mover events: {total_fm}')
print(f'Already had announcement date:    {already_had}')
print(f'Updated with researched dates:    {updated_count}')
print(f'Not found (no plant match):       {not_found}')
print(f'Events with announcement date:    {with_ann}/{total_fm} ({100*with_ann/total_fm:.0f}%)')
print(f'\nDate quality in lookup:')
print(f'  Firm (press release/filing/settlement): {firm_dates}')
print(f'  Approximate (month-level proxy):        {approx_dates}')
print(f'\nLead time analysis (announcement -> retirement):')
if leads:
    lead_days = [l[1] for l in leads]
    print(f'  Events with both dates: {len(leads)}')
    print(f'  Mean lead time:  {sum(lead_days)/len(lead_days):.0f} days '
          f'({sum(lead_days)/len(lead_days)/365:.1f} years)')
    print(f'  Median lead time: {sorted(lead_days)[len(lead_days)//2]:.0f} days')
    print(f'  Min: {min(lead_days)} days | Max: {max(lead_days)} days')
    # Show top 10 longest leads
    print(f'\n  Top 10 longest lead times:')
    for plant, days in sorted(leads, key=lambda x: -x[1])[:10]:
        print(f'    {plant:45s} {days:5d} days ({days/365:.1f} yr)')
    # Show any negative leads (announcement after retirement — data issue)
    negs = [(p, d) for p, d in leads if d < 0]
    if negs:
        print(f'\n  WARNING: {len(negs)} events with negative lead time:')
        for plant, days in negs:
            print(f'    {plant:45s} {days:5d} days')

print(f'\nSaved updated coal_retirement_events.csv')
