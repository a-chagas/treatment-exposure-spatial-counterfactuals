# Data used by the reproduction package

The reproduction scripts start from five processed inputs in
`data/processed/london_ulez_air/`:

1. `laqn_full_monthly_no2_2016_2024.csv` — monthly NO2 panel;
2. `laqn_station_ulez_assignment.csv` — monitoring-station coordinates, ULEZ membership,
   and station opening dates;
3. `ulez_2019.gpkg` — 2019 ULEZ boundary;
4. `ulez_2021.gpkg` — 2021 cumulative ULEZ boundary;
5. `ulez_2023.gpkg` — 2023 cumulative ULEZ boundary.

The R scripts reconstruct all exposure variables used in the paper from these inputs. In
particular, the static station-based exposure matrices are rebuilt on the outcome-monitor
universe, while the dynamic source-specific exposure matrices use the complete assignment
station universe, matching the empirical specifications reported in the manuscript.

## Source-data licensing and attribution

The London Air Quality Network (LAQN) data service operated by the Environmental Research
Group at Imperial College London makes its air-quality data available under the UK Open
Government Licence version 2.0 (OGL v2.0). The processed monitoring files in this package are
derived from those data.

Suggested attribution:

> Source: Environmental Research Group, Imperial College London, using data from the London
> Air Quality Network. Contains public sector information licensed under the Open Government
> Licence v2.0.

The 2019, 2021, and 2023 ULEZ boundary datasets distributed through the London Datastore are
also published under OGL v2.0 and identify Greater London Authority / Transport for London as
the relevant providers/authors.

Suggested attribution:

> Source: Greater London Authority / Transport for London, London Datastore. Contains public
> sector information licensed under the Open Government Licence v2.0.

OGL v2.0 permits copying, distribution, transmission, adaptation, and commercial or
non-commercial reuse subject to source attribution and the other licence conditions. The
licence is available at:

https://www.nationalarchives.gov.uk/doc/open-government-licence/version/2/

Source pages:

- LAQN API / data documentation: https://www.londonair.org.uk/Londonair/API/
- ULEZ 2019: https://data.london.gov.uk/dataset/ultra-low-emissions-zone-2019-v8onw
- ULEZ 2021: https://data.london.gov.uk/dataset/inner-ultra-low-emissions-zone-expansion-2021-vdjx4
- ULEZ 2023: https://data.london.gov.uk/dataset/london-wide-ultra-low-emission-zone-2023-vd455

These source-data terms are separate from the software licence for the reproduction code.
