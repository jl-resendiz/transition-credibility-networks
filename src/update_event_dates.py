"""Update coal_retirement_events.csv with exact dates from web research.

Dates sourced from utility press releases, news articles, EIA records,
GEM Wiki, Sierra Club, PJM deactivation notices, and company announcements.
Each date is tagged with confidence level and source.

Retirement dates are loaded from data/raw/events/manual_retirement_dates.csv,
a researcher-curated reference file with columns:
  plant_name, date, precision, source
"""
import csv, os
from _paths import raw_path, derived_path

# ── Load date lookup from CSV ─────────────────────────────────────────
# precision: 'day' = exact day known, 'month' = month known (use 15th)
# For plants with multiple units, use earliest unit deactivation date.

PLANT_DATES = {}
dates_csv = raw_path('events', 'manual_retirement_dates.csv')
with open(dates_csv, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        PLANT_DATES[row['plant_name']] = (
            row['date'], row['precision'], row['source']
        )


# ── Update the CSV ──────────────────────────────────────────────────
events_path = derived_path('events', 'coal_retirement_events.csv')

all_events = []
with open(events_path, 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    fieldnames = reader.fieldnames
    for row in reader:
        all_events.append(row)

updated_count = 0
already_had_date = 0
not_matched_fm = 0

for row in all_events:
    if row['is_first_mover'] != 'True' or row['is_matched'] != 'True':
        continue

    plant = row['plant_name']

    # Already has a date from wiki scraping
    if row.get('event_date') and row['event_date'].strip():
        already_had_date += 1
        continue

    # Look up in our research
    if plant in PLANT_DATES:
        date, precision, source = PLANT_DATES[plant]
        row['event_date'] = date
        updated_count += 1
    else:
        not_matched_fm += 1

# Write back
with open(events_path, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(all_events)

# Stats
total_fm = sum(1 for r in all_events if r['is_first_mover'] == 'True' and r['is_matched'] == 'True')
with_date = sum(1 for r in all_events if r['is_first_mover'] == 'True' and r['is_matched'] == 'True' and r.get('event_date') and r['event_date'].strip())
day_precision = sum(1 for plant, (d, p, s) in PLANT_DATES.items() if p == 'day')
month_precision = sum(1 for plant, (d, p, s) in PLANT_DATES.items() if p == 'month')

print(f'=== EVENT DATE UPDATE ===')
print(f'Total matched first-mover events: {total_fm}')
print(f'Already had date (wiki scraping): {already_had_date}')
print(f'Updated with web-researched dates: {updated_count}')
print(f'Still missing (plant name mismatch): {not_matched_fm}')
print(f'Events with date after update: {with_date}/{total_fm} ({100*with_date/total_fm:.0f}%)')
print(f'\nDate precision in lookup:')
print(f'  Day-level: {day_precision}')
print(f'  Month-level: {month_precision}')
print(f'\nSaved updated coal_retirement_events.csv')
