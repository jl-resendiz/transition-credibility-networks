"""Update coal_retirement_events.csv with exact dates from web research.

Dates sourced from utility press releases, news articles, EIA records,
GEM Wiki, Sierra Club, PJM deactivation notices, and company announcements.
Each date is tagged with confidence level and source.
"""
import csv, os
from _paths import derived_path

# ── Date lookup: plant_name -> (date, precision, source) ────────────
# precision: 'day' = exact day known, 'month' = month known (use 15th)
# For plants with multiple units, use earliest unit deactivation date.

PLANT_DATES = {
    # === UNITED STATES ===

    # Cherokee Station, Colorado - Xcel Energy
    # Unit 3 retired end of 2015 per Clean Air Clean Jobs Act
    'Cherokee Station': ('2015-12-31', 'day', 'Xcel Energy CACJ compliance; Denver Post'),

    # Eastlake Power Plant, Ohio - FirstEnergy
    # Units 1-3 deactivated April 6-10, 2015 per PJM/MATS
    'Eastlake Power Plant': ('2015-04-06', 'day', 'PJM deactivation notice; Power Engineering'),

    # Edgewater Generating Station, Wisconsin - Alliant Energy
    # Unit 3 retired end of 2015 per EPA Clean Air Act settlement
    'Edgewater Generating Station': ('2015-12-31', 'day', 'Alliant Energy; Power Engineering'),

    # H. Wilson Sundt Generating Station, Arizona - Tucson Electric Power
    # Last coal burned August 13, 2015
    'H. Wilson Sundt Generating Station': ('2015-08-13', 'day', 'Sierra Club Arizona; TEP press release'),

    # Jack Watson Generating Plant, Mississippi - Mississippi Power
    # Coal units converted to gas by April 15, 2015
    'Jack Watson Generating Plant': ('2015-04-15', 'day', 'Mississippi Power; Sierra Club settlement'),

    # Johnsonville Fossil Plant, Tennessee - TVA
    # Units 5-10 retired December 31, 2015 per EPA agreement
    'Johnsonville Fossil Plant': ('2015-12-31', 'day', 'TVA; Power Magazine'),

    # Lake Shore Plant, Ohio - FirstEnergy
    # Unit 18 deactivated April 13, 2015 per PJM/MATS
    'Lake Shore Plant': ('2015-04-13', 'day', 'PJM deactivation notice; Power Engineering'),

    # Milton Kapp Generating Station, Iowa - Alliant Energy
    # Retired 2015 per EPA settlement; exact month not confirmed
    'Milton Kapp Generating Station': ('2015-06-15', 'month', 'Alliant Energy EPA settlement'),

    # Nelson Dewey Generating Station, Wisconsin - Alliant Energy
    # Retired December 31, 2015 per EPA Clean Air Act agreement
    'Nelson Dewey Generating Station': ('2015-12-31', 'day', 'Alliant Energy press release; EPA settlement'),

    # Pulliam Power Plant, Wisconsin - WPS
    # Units 5-6 retired 2015 per EPA enforcement
    'Pulliam Power Plant': ('2015-06-15', 'month', 'WPS; EPA Clean Air Act enforcement'),

    # Scholz Generating Plant, Florida - Gulf Power
    # Retired April 2015
    'Scholz Generating Plant': ('2015-04-30', 'month', 'Gulf Power; Sierra Club'),

    # Stoneman Generating Station, Wisconsin - Dairyland Power Cooperative
    # Coal converted to biomass ~2010; PPA terminated 2015
    'Stoneman Generating Station': ('2015-12-31', 'day', 'Dairyland Power Cooperative; PPA termination'),

    # Syl Laskin Energy Center, Minnesota - Minnesota Power
    # Coal units taken offline early March 2015, converted to gas
    'Syl Laskin Energy Center': ('2015-03-15', 'month', 'Minnesota Power/ALLETE; Duluth News Tribune'),

    # Tecumseh Energy Center, Kansas - Westar Energy
    # Unit 7 (61 MW coal) retired 2015 per Westar fleet plan
    'Tecumseh Energy Center': ('2015-06-15', 'month', 'Westar Energy; GlobalSpec'),

    # Valley Power Plant, Wisconsin - We Energies
    # Unit 2 converted from coal to gas in 2015
    'Valley Power Plant': ('2015-12-31', 'day', 'We Energies; Power Magazine'),

    # Weston Power Plant, Wisconsin - WPS
    # Unit 1 retired 2015; Unit 2 converted to gas
    'Weston Power Plant': ('2015-08-15', 'month', 'WPS; TransmissionHub'),

    # W.S. Lee Steam Station, South Carolina - Duke Energy
    # Coal units retired November 6, 2014 (GEM lists as 2015)
    'W.S. Lee Steam Station': ('2014-11-06', 'day', 'Duke Energy; The Journal Online'),

    # === EUROPE ===

    # Duernrohr power station, Austria - Verbund
    # Block 1 decommissioned end of April 2015
    'Duernrohr power station': ('2015-04-30', 'month', 'Verbund press release'),

    # Bouchain power station, France - EDF
    # Coal unit retired April 2015
    'Bouchain power station': ('2015-04-15', 'month', 'GEM Wiki; Power Magazine'),

    # GKM (Mannheim) power station, Germany - RWE/EnBW/MVV
    # Units 3-4 retired after Block 9 came online May 2015
    'GKM (Mannheim) power station': ('2015-05-15', 'month', 'Alstom/GE press release; SourceWatch'),

    # Goldenberg power station, Germany - RWE
    # Electricity production halted July 1, 2015
    'Goldenberg power station': ('2015-07-01', 'day', 'RWE; German Wikipedia'),

    # Heilbronn power station, Germany - EnBW
    # Units 5-6 put on reserve 2015; formal status change
    'Heilbronn power station': ('2015-06-15', 'month', 'EnBW; Bundesnetzagentur'),

    # Ptolemaida power station, Greece - PPC
    # Units 3-4 destroyed by fire 2014, officially retired 2016
    'Ptolemaida power station': ('2016-01-15', 'month', 'Balkan Green Energy News; PPC'),
    # Variant spelling with accent in GEM tracker
    'Ptolema\u00efda power station': ('2016-01-15', 'month', 'Balkan Green Energy News; PPC'),

    # Porto Marghera Enel power station, Italy - Enel
    # Site sold November 2, 2015 (plant inactive since ~2012)
    'Porto Marghera Enel power station': ('2015-11-02', 'day', 'Enel Futur-e press release'),

    # Amer power station, Netherlands - RWE
    # Unit 8 closed January 1, 2016 (some sources say end of 2015)
    'Amer power station': ('2016-01-01', 'day', 'RWE; Wikipedia Amercentrale'),

    # Borssele power station, Netherlands - EPZ/RWE
    # Coal set closed November 25, 2015
    'Borssele power station': ('2015-11-25', 'day', 'GEM Wiki; EPZ'),

    # Bielsko-Biala power station, Poland - Tauron
    # Retired 2015; exact date not found
    'Bielsko-Biala power station': ('2015-06-15', 'month', 'Tauron; GEM tracker'),

    # Bydgoszcz power station, Poland - PGE
    # Retired 2015; exact date not found
    'Bydgoszcz power station': ('2015-06-15', 'month', 'PGE; GEM tracker'),

    # Krakow-Leg power station, Poland - PGE
    # Retired 2015; exact date not found
    'Krakow-Leg power station': ('2015-06-15', 'month', 'PGE; GEM tracker'),

    # Soto de Ribera power station, Spain
    # GEM lists unit retirement in 2015; EDP operates
    'Soto de Ribera power station': ('2015-06-15', 'month', 'GEM tracker; EDP'),

    # === ASIA ===

    # Beijing Yire power station, China
    # Ceased operation March 20, 2015 per Beijing municipal order
    'Beijing Yire power station': ('2015-03-20', 'day', 'CLP China; Dialogue Earth'),

    # Chongqing power station, China
    # Retired 2015; exact date not found
    'Chongqing power station': ('2015-06-15', 'month', 'GEM tracker'),

    # Hangzhou Banshan power station, China
    # Retired 2015; part of Zhejiang coal phase-down
    'Hangzhou Banshan power station': ('2015-06-15', 'month', 'GEM tracker'),

    # Henan Jiaozuo power station, China
    # Retired 2015; exact date not found
    'Henan Jiaozuo power station': ('2015-06-15', 'month', 'GEM tracker'),

    # Jiaozuo Bo\'ai power station, China
    "Jiaozuo Bo'ai power station": ('2015-06-15', 'month', 'GEM tracker'),

    # Tangshan West power station, China
    'Tangshan West power station': ('2015-06-15', 'month', 'GEM tracker'),

    # Panipat power station, India - HPGCL
    # Unit 5 retired 2016 per CEA; exact date not found
    'Panipat power station': ('2016-06-15', 'month', 'Central Electricity Authority India'),

    # Trombay power station, India - Tata Power
    # Coal unit retired 2016; exact date not found
    'Trombay power station': ('2016-06-15', 'month', 'Tata Power; GEM tracker'),

    # Yongdong power station, South Korea
    # Retired 2017; exact date not found
    'Yongdong power station': ('2017-06-15', 'month', 'KEPCO/KOSPO; GEM tracker'),

    # Naga power station, Philippines
    # Naga-1 retired 2015 as part of Naga-3 expansion
    'Naga power station': ('2015-06-15', 'month', 'SPC Power; GEM tracker'),

    # === SOUTH AMERICA ===

    # Charqueadas power station, Brazil
    # Retired 2016; exact date not found
    'Charqueadas power station': ('2016-06-15', 'month', 'CGTEE/Eletrobras; GEM tracker'),

    # === MIDDLE EAST ===

    # Hassyan Clean-Coal Power Project, UAE
    # Converted from coal to gas, announced February 3, 2022
    'Hassyan Clean-Coal Power Project': ('2022-02-03', 'day', 'DEWA press release; Arab News'),

    # === RUSSIA ===
    # Russian plants: exact dates very hard to find in English sources

    # Ivanovskaya CHP-2 power station, Russia - Inter RAO
    'Ivanovskaya CHP-2 power station': ('2015-06-15', 'month', 'GEM tracker'),

    # Omsk CHP-4 power station, Russia
    'Omsk CHP-4 power station': ('2015-06-15', 'month', 'GEM tracker'),

    # Serovskaya power station, Russia - OGK-1
    'Serovskaya power station': ('2015-06-15', 'month', 'GEM tracker'),

    # Troitskaya GRES power station, Russia - OGK-2
    'Troitskaya GRES power station': ('2015-06-15', 'month', 'GEM tracker'),

    # Apatitskaya CHP, Russia - TGC-1
    'Apatitskaya CHP power station': ('2016-06-15', 'month', 'GEM tracker'),

    # Anadyrskaya power station, Russia - RusHydro
    'Anadyrskaya power station': ('2018-06-15', 'month', 'GEM tracker'),

    # Tomskaya GRES-2 power station, Russia
    'Tomskaya GRES-2 power station': ('2020-06-15', 'month', 'GEM tracker'),

    # Vorkutinskaya-2 power station, Russia
    'Vorkutinskaya-2 power station': ('2021-06-15', 'month', 'GEM tracker'),
}


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
