--Understanding the data: total_cases/deaths is a daily tally of the total up to and including that day. True total is found on the last day of recording or can sum new cases/deaths columns. 
SELECT location, date, total_cases, new_cases, total_deaths, new_deaths
FROM "COVID_deaths"
WHERE location='Afghanistan'
ORDER BY date DESC;

--Looking at total cases vs total deaths over time
SELECT location, date, total_cases, total_deaths, ROUND(((CAST(total_deaths as decimal))/total_cases),2) AS death_percentage
FROM "COVID_deaths"
ORDER BY 1;

--Looking at perecentage deaths from total number of cases by using most recent date where total deaths and cases have been tallied up to that date.
SELECT location, date, total_cases, total_deaths, ROUND(((CAST(total_deaths as decimal))/total_cases),5)*100 AS death_percentage
FROM "COVID_deaths"
WHERE date = CAST('2022-08-28' as DATE);

--Looking at total deaths vs total cases by summing new cases and deaths columns instead of using most recent date.
SELECT location, SUM(new_cases) AS cases, SUM(new_deaths) AS deaths, ROUND((CAST (SUM(new_deaths)as decimal)/SUM(new_cases)),3)AS death_percentage
FROM "COVID_deaths"
GROUP BY location
ORDER BY 4 DESC;

--Getting just the value for the united states when not sure of exact location value; Percent of total cases that died
SELECT location, SUM(new_cases) AS cases, SUM(new_deaths) AS deaths, ROUND((CAST (SUM(new_deaths)as decimal)/SUM(new_cases)),5)*100 AS death_percentage
FROM "COVID_deaths"
WHERE location LIKE '% States'
GROUP BY location;

--Looking at the change in hopitalized and ICU patients over time and most current numbers
SELECT location,date, total_cases, hosp_patients, icu_patients, total_deaths
FROM "COVID_deaths"
WHERE location='United States';

---percentage of the population with confirmed cases for India, UK, South Korea, US, Australia and Sweden as of Aug 2022
SELECT location, date, total_cases,population, ROUND((CAST(total_cases as decimal)/population),5) * 100 AS percent_cases
FROM "COVID_deaths"
WHERE (location='Australia' OR location='United Kingdom' OR location='South Korea'OR location='India' OR location='United States' OR location='Sweden') AND date=CAST('08-28-2022' AS DATE)
GROUP BY 1,2,3,4;

--Percent of population that has died as of 8-28-2022
SELECT location, population,total_deaths,ROUND(((CAST(total_deaths as decimal))/population),5)*100 AS percentage_population_died
FROM "COVID_deaths"
WHERE date=(CAST('08-28-2022'AS DATE))
GROUP BY 1,2,3
ORDER BY 4 DESC;

--Percent of population that died for select countries
SELECT location, population,total_deaths,ROUND(((CAST(total_deaths as decimal))/population),5)*100 AS percentage_population_died
FROM "COVID_deaths"
WHERE (location= 'Canada'OR location='Australia' OR location='United Kingdom' OR location='South Korea'OR location='India' OR location='United States' OR location='Sweden') AND date=CAST('08-28-2022' AS DATE)
GROUP BY 1,2,3
ORDER BY 4 DESC;

---total global cases and deaths per continent
SELECT location AS continent, sum(new_cases) AS total_cases, sum(new_deaths) AS total_deaths
FROM "COVID_deaths"
where continent is null 
group by location;

---percent of global infections resulting in death per continent
SELECT location AS continent, ROUND(sum(cast(new_deaths as decimal))/sum(new_cases),4)*100 AS death_percentage
FROM "COVID_deaths"
where continent is null and (location!='High income' AND location!='Low income' AND location!='Lower middle income' AND location!='Upper middle income' AND location!='International' AND location!='World')
group by location;

---percent of continent population that died from COVID
select location, population, ROUND(sum(CAST(new_deaths as decimal))/population,5)*100 AS death_percent
from "COVID_deaths"
where continent is null and (location!='High income' AND location!='Low income' AND location!='Lower middle income' AND location!='Upper middle income' AND location!='International' AND location!='World')
group by location, population
ORDER BY 3 DESC;

---Checking out data in vaccines table
SELECT *
From "COVID_vax"
LIMIt 10;

---Evaluating the effect of vaccines on case numbers
SELECT "COVID_deaths".location,"COVID_deaths".date,"COVID_deaths".total_cases,"COVID_vax".people_vaccinated,"COVID_vax".people_fully_vaccinated
FROM "COVID_deaths"
JOIN "COVID_vax"
  ON "COVID_vax"."location" = "COVID_deaths"."location" and
  "COVID_vax".date = "COVID_deaths".date
ORDER BY 1,2;

---Does the number of new deaths decline after large portion of population is fully vaccinated?
SELECT  "COVID_deaths".location,"COVID_deaths".date, ROUND(((CAST("COVID_vax".people_fully_vaccinated AS DECIMAL)/"COVID_vax".population)*100),4) AS percent_vaccinated, "COVID_deaths".new_deaths 
FROM "COVID_deaths"
JOIN "COVID_vax"
  ON "COVID_vax"."location" = "COVID_deaths"."location" and
  "COVID_vax".date = "COVID_deaths".date
WHERE ("COVID_deaths".date BETWEEN (CAST('2021-01-01' as DATE)) AND (CAST('2021-06-01' as DATE)) OR "COVID_deaths".date BETWEEN (CAST('2022-01-01' AS DATE)) AND (CAST('2022-06-01' AS DATE))) AND ("COVID_deaths".location='United States' OR "COVID_deaths".location='United Kingdom'OR "COVID_deaths".location='Sweden' OR "COVID_deaths".location='South Korea' OR "COVID_deaths".location='Australia')
ORDER BY 1,2

---Pattern of case numbers and deaths by level of diabetes prevalence(percent of pop with diabetes)
SELECT "COVID_deaths".location,"COVID_deaths".date,ROUND((CAST("COVID_deaths".total_deaths as decimal)/"COVID_deaths".total_cases),5)*100 AS death_percentage,"COVID_vax".diabetes_prevalence
FROM "COVID_deaths"
JOIN "COVID_vax"
  ON "COVID_vax"."location" = "COVID_deaths"."location" and
  "COVID_vax".date = "COVID_deaths".date
WHERE "COVID_deaths".date=CAST('08-28-2022' AS DATE) AND "COVID_vax".diabetes_prevalence IS NOT NULL
ORDER BY 3;

---Do locations with more smokers have increased infections resulting in death?
SELECT "COVID_deaths".location, "COVID_deaths".total_cases,("COVID_vax".female_smokers + "COVID_vax".male_smokers) AS number_smokers,ROUND((CAST("COVID_deaths".total_deaths as decimal)/"COVID_deaths".total_cases),5)*100 AS death_percentage
FROM "COVID_deaths"
JOIN "COVID_vax"
  ON "COVID_vax"."location" = "COVID_deaths"."location" and
  "COVID_vax".date = "COVID_deaths".date
WHERE "COVID_deaths".date=CAST('08-28-2022' AS DATE)
GROUP BY 1,2,3,4
ORDER BY 4 DESC;

---Did policy stringency effect COVID death percentage in population? Stringency index is 1-100 with 100 being strictest.
--Find max and min stringency for each location
SELECT location, MAX(stringency_index),MIN(stringency_index)
FROM "COVID_vax"
GROUP BY location;
--Found that each location does have a MAX and MIN so will take the average stringecy to see if there is much variation
SELECT location, AVG(stringency_index)
FROM "COVID_vax"
GROUP BY location
ORDER BY 2;
---Didn't find much variation in average stringency so changed plan to get the change in stringency over time and the change in percent of cases resulting in death over time for select countries known to have strict vs. loose stringency approaches
SELECT "COVID_deaths".date,"COVID_deaths".location,ROUND((CAST("COVID_deaths".total_deaths as decimal))/"COVID_deaths".total_cases,5) AS percent_deaths, "COVID_vax".stringency_index,
  CASE
    WHEN "COVID_vax".stringency_index BETWEEN 1.0 AND 25.999 THEN 'Low'
	WHEN "COVID_vax".stringency_index BETWEEN 26.0 AND 75.999 THEN 'Moderate'
	WHEN "COVID_vax".stringency_index BETWEEN 76.0 AND 100.0 THEN 'HIGH'
	ELSE 'None'
  END Stringency_Measures
FROM "COVID_deaths"
JOIN "COVID_vax"
ON "COVID_vax".continent = "COVID_deaths".continent and "COVID_vax".date = "COVID_deaths".date
WHERE "COVID_deaths".location='United States' OR "COVID_deaths".location='Australia' OR "COVID_deaths".location='Canada' OR "COVID_deaths".location='United Kingdom' OR "COVID_deaths".location='Sweden'OR "COVID_deaths".location='South Korea'
GROUP BY 1,2,3,4
ORDER BY "COVID_deaths".location;

----How does viral reproduction rate change over time in the US?
SELECT date, location, reproduction_rate
FROM "COVID_deaths"
WHERE location='United States';

--Does viral reproduction rate change following increased or decreased stringency measures?
SELECT "COVID_deaths".date,"COVID_deaths".location,"COVID_deaths".reproduction_rate, "COVID_vax".stringency_index,
  CASE
    WHEN "COVID_vax".stringency_index BETWEEN 1.0 AND 25.999 THEN 'Low'
	WHEN "COVID_vax".stringency_index BETWEEN 26.0 AND 75.999 THEN 'Moderate'
	WHEN "COVID_vax".stringency_index BETWEEN 76.0 AND 100.0 THEN 'HIGH'
	ELSE 'None'
  END Stringency_Measures
FROM "COVID_deaths"
JOIN "COVID_vax"
ON "COVID_vax".continent = "COVID_deaths".continent and "COVID_vax".date = "COVID_deaths".date
WHERE "COVID_deaths".location='United States' OR "COVID_deaths".location='Australia' OR "COVID_deaths".location='Canada' OR "COVID_deaths".location='United Kingdom' OR "COVID_deaths".location='Sweden'OR "COVID_deaths".location='South Korea'
GROUP BY 1,2,3,4
ORDER BY "COVID_deaths".location;
---Creating some views for later visualization
CREATE VIEW COVID_law_strigency_and_percent_of_cases_resulting_in_deaths_over_time as
SELECT "COVID_deaths".date,"COVID_deaths".location,ROUND((CAST("COVID_deaths".total_deaths as decimal))/"COVID_deaths".total_cases,5) AS percent_deaths, "COVID_vax".stringency_index,
  CASE
    WHEN "COVID_vax".stringency_index BETWEEN 1.0 AND 25.999 THEN 'Low'
	WHEN "COVID_vax".stringency_index BETWEEN 26.0 AND 75.999 THEN 'Moderate'
	WHEN "COVID_vax".stringency_index BETWEEN 76.0 AND 100.0 THEN 'HIGH'
	ELSE 'None'
  END Stringency_Measures
FROM "COVID_deaths"
JOIN "COVID_vax"
ON "COVID_vax".continent = "COVID_deaths".continent and "COVID_vax".date = "COVID_deaths".date
WHERE "COVID_deaths".location='United States' OR "COVID_deaths".location='Australia' OR "COVID_deaths".location='Canada' OR "COVID_deaths".location='United Kingdom' OR "COVID_deaths".location='Sweden'OR "COVID_deaths".location='South Korea'
GROUP BY 1,2,3,4
ORDER BY "COVID_deaths".location;

--View 2
CREATE VIEW covid_law_stringency_viral_reproduction_rate AS
SELECT "COVID_deaths".date,"COVID_deaths".location,"COVID_deaths".reproduction_rate, "COVID_vax".stringency_index,
  CASE
    WHEN "COVID_vax".stringency_index BETWEEN 1.0 AND 25.999 THEN 'Low'
	WHEN "COVID_vax".stringency_index BETWEEN 26.0 AND 75.999 THEN 'Moderate'
	WHEN "COVID_vax".stringency_index BETWEEN 76.0 AND 100.0 THEN 'HIGH'
	ELSE 'None'
  END Stringency_Measures
FROM "COVID_deaths"
JOIN "COVID_vax"
ON "COVID_vax".continent = "COVID_deaths".continent and "COVID_vax".date = "COVID_deaths".date
WHERE "COVID_deaths".location='United States' OR "COVID_deaths".location='Australia' OR "COVID_deaths".location='Canada' OR "COVID_deaths".location='United Kingdom' OR "COVID_deaths".location='Sweden'OR "COVID_deaths".location='South Korea'
GROUP BY 1,2,3,4
ORDER BY "COVID_deaths".location;

--View 3 
CREATE VIEW diabetes_prevalence_and_death_rate AS
SELECT "COVID_deaths".location,"COVID_deaths".date,ROUND((CAST("COVID_deaths".total_deaths as decimal)/"COVID_deaths".total_cases),5)*100 AS death_percentage,"COVID_vax".diabetes_prevalence
FROM "COVID_deaths"
JOIN "COVID_vax"
  ON "COVID_vax"."location" = "COVID_deaths"."location" and
  "COVID_vax".date = "COVID_deaths".date
WHERE "COVID_deaths".date=CAST('08-28-2022' AS DATE) AND "COVID_vax".diabetes_prevalence IS NOT NULL
ORDER BY 3;

--View 4
CREATE VIEW change_in_case_numbers_as_vaccinations_increased AS
SELECT "COVID_deaths".location,"COVID_deaths".date,"COVID_deaths".total_cases,"COVID_vax".people_vaccinated,"COVID_vax".people_fully_vaccinated
FROM "COVID_deaths"
JOIN "COVID_vax"
  ON "COVID_vax"."location" = "COVID_deaths"."location" and
  "COVID_vax".date = "COVID_deaths".date
ORDER BY 1,2;

--View 5
CREATE VIEW number_of_new_deaths_with_increasing_vaccination AS
SELECT  "COVID_deaths".location,"COVID_deaths".date, ROUND(((CAST("COVID_vax".people_fully_vaccinated AS DECIMAL)/"COVID_vax".population)*100),4) AS percent_vaccinated, "COVID_deaths".new_deaths 
FROM "COVID_deaths"
JOIN "COVID_vax"
  ON "COVID_vax"."location" = "COVID_deaths"."location" and
  "COVID_vax".date = "COVID_deaths".date
WHERE ("COVID_deaths".date BETWEEN (CAST('2021-01-01' as DATE)) AND (CAST('2021-06-01' as DATE)) OR "COVID_deaths".date BETWEEN (CAST('2022-01-01' AS DATE)) AND (CAST('2022-06-01' AS DATE))) AND ("COVID_deaths".location='United States' OR "COVID_deaths".location='United Kingdom'OR "COVID_deaths".location='Sweden' OR "COVID_deaths".location='South Korea' OR "COVID_deaths".location='Australia')
ORDER BY 1,2;

--View 6
CREATE VIEW percentage_of_deaths_per_country AS
SELECT location, population,total_deaths,ROUND(((CAST(total_deaths as decimal))/population),5)*100 AS percentage_population_died
FROM "COVID_deaths"
WHERE (location= 'Canada'OR location='Australia' OR location='United Kingdom' OR location='South Korea'OR location='India' OR location='United States' OR location='Sweden') AND date=CAST('08-28-2022' AS DATE)
GROUP BY 1,2,3
ORDER BY 4 DESC;

--View 7
CREATE VIEW total_cases_and_deaths_per_continent AS
SELECT location AS continent, sum(new_cases) AS total_cases, sum(new_deaths) AS total_deaths
FROM "COVID_deaths"
where continent is null 
group by location;