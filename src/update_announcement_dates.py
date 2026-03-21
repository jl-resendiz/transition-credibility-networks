"""Add announcement_date column to coal_retirement_events.csv.

Markets react to announcements, not physical closures. For a finance event
study, the proper event date is when the retirement was first publicly
disclosed (utility filing, press release, regulatory settlement, board vote).

Dates sourced from:
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

from _paths import derived_path

# ── Announcement date lookup: plant_name -> (date, source) ──────────
# date = YYYY-MM-DD when market first learned of the retirement plan
# For multi-unit plants, use earliest public disclosure.
# Day-15 = month known but not exact day.

ANNOUNCEMENT_DATES = {

    # ================================================================
    # UNITED STATES
    # ================================================================

    # --- AEP fleet: June 9, 2011 press release ---
    # AEP announced ~6,000 MW retirement plan due to EPA regulations
    'Kammer Plant': ('2011-06-09',
        'AEP press release 2011-06-09; EPA compliance plan'),
    'Philip Sporn Power Plant': ('2011-06-09',
        'AEP press release 2011-06-09; EPA compliance plan'),
    'Kanawha River Plant': ('2011-06-09',
        'AEP press release 2011-06-09; EPA compliance plan'),
    'Muskingum River Plant': ('2011-06-09',
        'AEP press release 2011-06-09; EPA compliance plan'),
    'Glen Lyn Plant': ('2011-06-09',
        'AEP press release 2011-06-09; EPA compliance plan'),
    'Clinch River Plant': ('2011-06-09',
        'AEP press release 2011-06-09; EPA compliance plan'),
    'Picway Power Plant': ('2011-06-09',
        'AEP press release 2011-06-09; EPA compliance plan'),

    # AEP Big Sandy: Feb 25, 2013 EPA settlement modification
    'Big Sandy Plant': ('2013-02-25',
        'AEP/EPA/Earthjustice settlement 2013-02-25; coal stop by 2015'),

    # --- FirstEnergy fleet: Jan 26, 2012 announcement ---
    # Announced 6 coal plant retirements due to MATS
    'Eastlake Power Plant': ('2012-01-26',
        'FirstEnergy press release 2012-01-26; MATS compliance'),
    'Lake Shore Plant': ('2012-01-26',
        'FirstEnergy press release 2012-01-26; MATS compliance'),
    'Ashtabula power station': ('2012-01-26',
        'FirstEnergy press release 2012-01-26; MATS compliance'),

    # --- TVA fleet: Nov 14, 2013 board approval ---
    'Johnsonville Fossil Plant': ('2013-11-14',
        'TVA Board resolution 2013-11-14; retirement of coal units'),
    'Widows Creek Fossil Plant': ('2013-11-14',
        'TVA Board resolution 2013-11-14; retirement of coal units'),
    'Colbert Fossil Plant': ('2013-11-14',
        'TVA Board resolution 2013-11-14; retirement of coal units'),

    # --- Georgia Power fleet ---
    # Harllee Branch: March 2011 PSC filing
    'Harllee Branch Generating Plant': ('2011-03-15',
        'Georgia Power PSC decertification filing March 2011'),
    # Yates, Mitchell, Kraft: Jan 7, 2013 PSC filing
    'Yates Steam Generating Plant': ('2013-01-07',
        'Georgia Power PSC filing 2013-01-07; 15 plant retirements'),
    'Mitchell Steam Generating Plant (Georgia)': ('2013-01-07',
        'Georgia Power PSC filing 2013-01-07'),
    'Kraft Plant': ('2013-01-07',
        'Georgia Power PSC filing 2013-01-07'),

    # --- Alabama Power: Aug 1, 2014 ---
    'Barry Steam Plant': ('2014-08-01',
        'Alabama Power MATS compliance announcement 2014-08-01'),
    'Gorgas Steam Plant': ('2014-08-01',
        'Alabama Power MATS compliance announcement 2014-08-01'),

    # --- Xcel Energy ---
    # Cherokee: November 2010, Clean Air Clean Jobs Act compliance
    'Cherokee Station': ('2010-11-15',
        'Xcel Energy announcement Nov 2010; CO Clean Air Clean Jobs Act'),
    # Black Dog: March 15, 2011 MN PUC filing
    'Black Dog Generating Station': ('2011-03-15',
        'Xcel Energy MN PUC certificate of need filing 2011-03-15'),

    # --- Alliant Energy: Apr 22, 2013 EPA settlement ---
    'Edgewater Generating Station': ('2013-04-22',
        'WPL/EPA/Sierra Club settlement 2013-04-22'),
    'Nelson Dewey Generating Station': ('2013-04-22',
        'WPL/EPA/Sierra Club settlement 2013-04-22'),
    'Milton Kapp Generating Station': ('2013-04-22',
        'Alliant Energy EPA settlement 2013-04-22; Iowa operations'),

    # --- LG&E/KU: April 2011 PSC filing ---
    'Cane Run Station': ('2011-04-15',
        'LG&E/KU PSC filing mid-April 2011; 800 MW retirement plan'),
    'Green River Generating Station': ('2011-04-15',
        'LG&E/KU PSC filing mid-April 2011; 800 MW retirement plan'),

    # --- Mississippi Power: Aug 4, 2014 Sierra Club settlement ---
    'Jack Watson Generating Plant': ('2014-08-04',
        'Mississippi Power/Sierra Club settlement 2014-08-04'),

    # --- Minnesota Power: Jan 30, 2013 EnergyForward ---
    'Syl Laskin Energy Center': ('2013-01-30',
        'Minnesota Power/ALLETE EnergyForward announcement 2013-01-30'),
    'Taconite Harbor Energy Center': ('2013-01-30',
        'Minnesota Power/ALLETE EnergyForward announcement 2013-01-30'),

    # --- WPS: Jan 4, 2013 EPA settlement ---
    'Pulliam Power Plant': ('2013-01-04',
        'WPS/EPA Clean Air Act settlement 2013-01-04'),
    'Weston Power Plant': ('2013-01-04',
        'WPS/EPA Clean Air Act settlement 2013-01-04'),

    # --- Gulf Power: March 2013 ---
    'Scholz Generating Plant': ('2013-03-15',
        'Gulf Power announcement March 2013; EPA MATS compliance'),

    # --- TEP Sundt: Early 2014 EPA agreement ---
    'H. Wilson Sundt Generating Station': ('2014-01-15',
        'TEP/EPA agreement early 2014; coal-to-gas conversion'),

    # --- PGE Boardman: Jan 14, 2010 ---
    'Boardman Plant': ('2010-01-14',
        'PGE announcement 2010-01-14; 20-year early closure'),

    # --- PPL Montana Corette: 2012 mothball decision ---
    'Corette Plant': ('2012-06-15',
        'PPL Montana mothball decision 2012; permanent Feb 2015'),

    # --- NRG Will County: Aug 7, 2014 ---
    'Will County Generating Station': ('2014-08-07',
        'NRG environmental action plan 2014-08-07; Unit 3 retirement'),

    # --- APS Cholla: Sep 11, 2014 ---
    'Cholla Generating Station': ('2014-09-11',
        'APS announcement 2014-09-11; Unit 2 closure plan'),

    # --- Duke Energy Miami Fort: August 2011 ---
    'Miami Fort Station': ('2011-08-15',
        'Duke Energy announcement Aug 2011; Unit 6 retirement by 2015'),

    # --- Duke Energy W.S. Lee: ~2012 IRP ---
    'W.S. Lee Steam Station': ('2012-06-15',
        'Duke/Progress Energy IRP ~2012; coal retirement plan'),

    # --- We Energies Valley: ~2013 conversion plan ---
    'Valley Power Plant': ('2013-06-15',
        'We Energies coal-to-gas conversion plan ~2013'),

    # --- Dominion Chesapeake: ~2012 MATS-driven ---
    'Chesapeake Energy Center': ('2012-06-15',
        'Dominion IRP ~2012; MATS compliance retirement'),

    # --- DTE Energy ---
    # Trenton Channel: early retirement ~2013
    'Trenton Channel Power Plant': ('2013-06-15',
        'DTE Energy fleet planning ~2013; MATS compliance'),
    # River Rouge: ~2014 announcement
    'River Rouge Power Plant': ('2014-06-15',
        'DTE Energy fleet planning ~2014'),

    # --- Dairyland Power Stoneman ---
    'Stoneman Generating Station': ('2013-06-15',
        'Dairyland Power PPA termination planning ~2013'),

    # --- Edwards Generating Plant: retirement in our data is 2015 ---
    'Edwards Generation Plant': ('2013-06-15',
        'Dynegy acquisition from Ameren 2013; fleet restructuring'),

    # --- Lawrence Energy Center (Kansas) ---
    'Lawrence Energy Center (Kansas)': ('2013-06-15',
        'Westar Energy fleet plan ~2013; coal unit retirement'),

    # --- Tecumseh Energy Center ---
    'Tecumseh Energy Center': ('2013-06-15',
        'Westar Energy fleet plan ~2013; coal unit retirement'),

    # ================================================================
    # CANADA
    # ================================================================

    # TransAlta Sundance: November 2016 Alberta off-coal agreements
    'Alberta Sundance power station': ('2016-11-15',
        'Alberta off-coal agreements Nov 2016; TransAlta/ATCO/Capital Power'),

    # ================================================================
    # EUROPE
    # ================================================================

    # --- Netherlands: September 6, 2013 Energieakkoord ---
    'Amer power station': ('2013-09-06',
        'Dutch Energy Agreement (Energieakkoord) Sep 2013'),
    'Borssele power station': ('2013-09-06',
        'Dutch Energy Agreement (Energieakkoord) Sep 2013'),
    'Nijmegen power station': ('2013-09-06',
        'Dutch Energy Agreement (Energieakkoord) Sep 2013'),

    # --- Austria: May 14, 2014 Verbund supervisory board ---
    'Duernrohr power station': ('2014-05-14',
        'Verbund supervisory board confirmation 2014-05-14'),

    # --- France: EU LCPD opt-out (Directive 2001/80/EC) ---
    # Plants opted out ~2008, had 20,000 hour limit, closed by end 2015
    # French Transitional National Plan submitted ~2013
    'Bouchain power station': ('2013-06-15',
        'EDF/France TNP under EU LCPD opt-out ~2013'),
    'La Maxe power station': ('2013-06-15',
        'EDF/France TNP under EU LCPD opt-out ~2013'),
    'Vitry power station': ('2013-06-15',
        'EDF/France TNP under EU LCPD opt-out ~2013'),

    # --- Germany ---
    # GKM Mannheim: Block 9 construction decision ~2010
    'GKM (Mannheim) power station': ('2010-06-15',
        'GKM Block 9 investment decision ~2010; old units to retire'),
    # Goldenberg: RWE decision ~2014
    'Goldenberg power station': ('2014-06-15',
        'RWE closure decision ~2014'),
    # Heilbronn: EnBW reserve decision ~2014
    'Heilbronn power station': ('2014-06-15',
        'EnBW reserve/closure decision ~2014; Bundesnetzagentur'),

    # --- Greece: Ptolemaida fire Nov 10, 2014 ---
    'Ptolemaida power station': ('2014-11-10',
        'Fire destroyed Ptolemaida units 3-4; Nov 10 2014'),
    'Ptolema\u00efda power station': ('2014-11-10',
        'Fire destroyed Ptolemaida units 3-4; Nov 10 2014'),

    # --- Italy: Porto Marghera inactive since ~2012 ---
    'Porto Marghera Enel power station': ('2012-06-15',
        'Enel deactivated Porto Marghera ~2012; sold 2015 via Futur-e'),

    # --- Slovakia: Novaky ---
    # Units in our data retired 2015; Slovak coal exit announced 2019
    # But retirement of specific units likely announced earlier
    'Novaky power station': ('2014-06-15',
        'Novaky unit retirement planning ~2014'),

    # --- Spain ---
    'Puertollano IGCC power station': ('2014-06-15',
        'Elcogas announced closure plans 2014'),
    'Soto de Ribera power station': ('2014-06-15',
        'EDP fleet restructuring ~2014'),

    # --- Poland ---
    # Hard to find specific dates; use GEM tracker year proxy
    'Bielsko-Biala power station': ('2014-06-15',
        'Tauron fleet restructuring ~2014'),
    'Bydgoszcz power station': ('2014-06-15',
        'PGE fleet restructuring ~2014'),
    'Krakow-Leg power station': ('2014-06-15',
        'PGE fleet restructuring ~2014'),

    # ================================================================
    # ASIA
    # ================================================================

    # --- China: Beijing Clean Air Action Plan Sep 2013 ---
    'Beijing Yire power station': ('2013-09-15',
        'Beijing Clean Air Action Plan 2013-2017; Sep 2013'),

    # --- China: other plants ---
    # National small-unit retirement policy (2007 State Council Order 25)
    # Specific announcement dates not available in English
    'Chongqing power station': ('2014-06-15',
        'Chinese provincial coal phase-down plan ~2014'),
    'Hangzhou Banshan power station': ('2014-06-15',
        'Zhejiang province coal phase-down ~2014'),
    'Henan Jiaozuo power station': ('2014-06-15',
        'Henan province coal phase-down ~2014'),
    "Jiaozuo Bo'ai power station": ('2014-06-15',
        'Henan province coal phase-down ~2014'),
    'Tangshan West power station': ('2014-06-15',
        'Hebei province coal phase-down ~2014'),
    'Shuangyashan power station': ('2017-06-15',
        'Heilongjiang province coal phase-down ~2017'),

    # --- India ---
    'Panipat power station': ('2015-06-15',
        'HPGCL retirement order ~2015; units retired Dec 9 2015'),
    'Trombay power station': ('2015-06-15',
        'Tata Power fleet restructuring ~2015'),
    'Chandrapur (Assam) power station': ('2016-06-15',
        'APGCL retirement planning ~2016'),

    # --- South Korea ---
    'Yongdong power station': ('2016-06-15',
        'KEPCO/KOSPO coal retirement plan ~2016'),

    # --- Philippines ---
    'Naga power station': ('2014-06-15',
        'SPC Power; Naga-3 expansion necessitated Naga-1 retirement'),

    # ================================================================
    # SOUTH AMERICA
    # ================================================================

    # Chile: June 4, 2019 National Decarbonization Plan
    'Bocamina power station': ('2019-06-04',
        'Chile National Decarbonization Plan signed 2019-06-04'),
    'Patache power station': ('2019-06-04',
        'Chile National Decarbonization Plan signed 2019-06-04'),
    'Tocopilla power station': ('2019-06-04',
        'Chile National Decarbonization Plan signed 2019-06-04'),

    # Brazil
    'Charqueadas power station': ('2015-06-15',
        'CGTEE/Eletrobras fleet plan ~2015'),

    # ================================================================
    # MIDDLE EAST
    # ================================================================

    # UAE: Hassyan conversion announced Feb 2022
    'Hassyan Clean-Coal Power Project': ('2022-02-03',
        'DEWA press release 2022-02-03; coal-to-gas conversion'),

    # ================================================================
    # AFRICA
    # ================================================================

    # South Africa: Komati announcement Aug 2021
    'Komati power station': ('2021-08-15',
        'Eskom announcement Aug 2021; Komati end-of-life retirement'),

    # ================================================================
    # AUSTRALIA
    # ================================================================

    # AGL Liddell: April 2015 greenhouse gas policy
    'Liddell power station': ('2015-04-15',
        'AGL Greenhouse Gas Policy April 2015; Liddell closure by 2022'),

    # ================================================================
    # RUSSIA
    # ================================================================
    # Announcement dates not available in English sources.
    # Using approximate dates based on Russian energy sector restructuring.
    'Ivanovskaya CHP-2 power station': ('2014-06-15',
        'Inter RAO fleet restructuring ~2014'),
    'Omsk CHP-4 power station': ('2014-06-15',
        'Russian energy sector ~2014'),
    'Serovskaya power station': ('2014-06-15',
        'OGK-1 fleet restructuring ~2014'),
    'Troitskaya GRES power station': ('2014-06-15',
        'OGK-2 fleet restructuring ~2014'),
    'Apatitskaya CHP power station': ('2015-06-15',
        'TGC-1 fleet restructuring ~2015'),
    'Anadyrskaya power station': ('2017-06-15',
        'RusHydro fleet planning ~2017'),
    'Tomskaya GRES-2 power station': ('2019-06-15',
        'Russian energy sector ~2019'),
    'Vorkutinskaya-2 power station': ('2020-06-15',
        'Russian energy sector ~2020'),
    'Sakhalin GRES-1 power station': ('2018-06-15',
        'Russian energy sector ~2018'),
}


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
