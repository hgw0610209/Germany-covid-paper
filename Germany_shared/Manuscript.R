#################################################################################
# This Rscript contains:
# (1) Running the health models in paper to obtain Table 3.
# (2) Running the health models in paper to obtain Table 4.
# (3) Fit model to up-to-date COVID-19 data.
#################################################################################
library(spdep)
library(INLA) 
library(dplyr)
library(sf)
library(stringr)

# load in the shapefile with full dataset
load(file="Final_shp.RData")

# create neighbourhood matrix
W <- poly2nb(Final_shp) 
W.moran <- nb2listw(W,zero.policy = TRUE)
W <- nb2mat(W, style="B",zero.policy = TRUE)

# the id for random effects
model_data <- Final_shp@data
model_data$ID <- 1:nrow(model_data)


##########################################(1) Running the health models in paper to obtain Table 3.###################################

### fit the health models

spatialModel_SIR <- vector(mode = "list", length = 4)

# set priors
prior <- list(
    theta1 = list(
        prior = "pc.prec",
        param = c(1, 0.5)),
    theta2 = list(
        prior = "gaussian",
        param = c(0, 1.8))
)

 set.seed(365)
# leroux model
spatialModel_SIR[[1]] <- inla(cumul_cases ~ no2+pm25+so2+Temperature+scale(Benzene)+scale(Aresenic)+scale(Cadmium)+scale(Nickel)+scale(popDensity)+f(ID, model = 'besagproper2', graph = W, hyper = prior), E = E, family = 'poisson', data = model_data
                    ,control.compute = list(config=TRUE,dic = TRUE, waic = TRUE)
                    ,control.predictor = list(compute = TRUE))
# bym
spatialModel_SIR[[2]] <- inla(cumul_cases ~ no2+pm25+so2+Temperature+scale(Benzene)+scale(Aresenic)+scale(Cadmium)+scale(Nickel)+scale(popDensity)+f(ID, model = 'bym', graph = W), E = E, family = 'poisson', data = model_data
                    ,control.compute = list(config=TRUE,dic = TRUE, waic = TRUE)
                    ,control.inla = list(h=0.00001)
                    ,control.predictor = list(compute = TRUE))
# besag
spatialModel_SIR[[3]] <- inla(cumul_cases ~ no2+pm25+so2+Temperature+scale(Benzene)+scale(Aresenic)+scale(Cadmium)+scale(Nickel)+scale(popDensity)+f(ID, model = 'besag', graph = W), E = E, family = 'poisson', data = model_data
                    ,control.compute = list(config=TRUE,dic = TRUE, waic = TRUE)
                    ,control.predictor = list(compute = TRUE))
# iid
spatialModel_SIR[[4]] <- inla(cumul_cases ~ no2+pm25+so2+Temperature+scale(Benzene)+scale(Aresenic)+scale(Cadmium)+scale(Nickel)+scale(popDensity)+f(ID, model = 'iid', graph = W), E = E, family = 'poisson', data = model_data
                    ,control.compute = list(config=TRUE,dic = TRUE, waic = TRUE)
                    ,control.predictor = list(compute = TRUE))

performanceTable <- mapply(function(x){round(x$waic$waic,2)}, spatialModel_SIR, SIMPLIFY = T)

resultTable <- mapply(function(x){round(cbind(100*(exp(x$summary.fixed[-1,c(grep("0.5quant", names(x$summary.fixed)),grep("0.025|0.975", names(x$summary.fixed))),])-1),
                                              1- mapply(function(xx){inla.pmarginal(0,xx)}, x$marginals.fixed, SIMPLIFY =T)[-1]
),2)[,-1]}, spatialModel_SIR, SIMPLIFY = F)

resultTable <- mapply(function(x){
  characterTable <- NULL
for(i in 1:(ncol(x)-1))
{
    characterTable <- cbind(characterTable,str_pad(sprintf("%.2f",round(x[,i],2)),6,pad=" "))
}
   characterTable <- cbind(characterTable,str_pad(sprintf("%.2f",round(x[,ncol(x)],2)),3,pad=" "))
   
 cbind(characterTable[,1], paste0("(",characterTable[,2],", ",characterTable[,3],")"), paste0("[",characterTable[,4],"]"))
  
},resultTable)

resultTable <- matrix(resultTable, nrow = 9)

resultTable <- rbind(resultTable, as.character(c("", performanceTable[1],"","",performanceTable[2],"","",performanceTable[3],"","",performanceTable[4],"")))

rownames(resultTable) <- NULL
rownames(resultTable) <- c("NO2","PM25","SO2","Temperature","Benzene","Aresenic","Cadmium","Nickel", "popDensity","WAIC")

# Table 3 in manuscript
resultTable

##########################################(2) Running the health models in paper to obtain Table 4.###################################

# Just rerun the above code by replacing "pm25" by "pm10" in the formula, will get the result for Table 4. 

##########################################(3) Fit model to up-to-date COVID-19 data.#################################################
## note that parts of the R code are original from kaggle (author: Heads or Tails), which has been cited in the paper.  
## https://www.kaggle.com/headsortails/covid19-tracking-germany
library(rstan)
library(robustbase)
library(data.table)
library(lubridate)
library(dplyr)
library(sp)
library(rgdal)
library(here)
library(INLA)
library(Epi)
library(survival)
library(denstrip)
library(MASS) 

#https://www.kaggle.com/headsortails/covid19-tracking-germany
libs <- c('dplyr', 'tibble',      # wrangling
          'stringr', 'readr',     # strings, input
          'lubridate', 'tidyr',   # time, wrangling
          'knitr', 'kableExtra',  # table styling
          'ggplot2', 'viridis',   # visuals
          'gganimate', 'sf',      # animations, maps
          'ggthemes')             # visuals
invisible(lapply(libs, library, character.only = TRUE))

infile <- "https://opendata.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0.csv"
covid_de <- read_csv(infile, col_types = cols())


covid_de <- covid_de %>% 
  dplyr::select(state = Bundesland,
                county = Landkreis,
                age_group = Altersgruppe,
                gender = Geschlecht,
                cases = AnzahlFall,
                deaths = AnzahlTodesfall,
                recovered = AnzahlGenesen,
                date = Meldedatum) %>% 
  mutate(date = date(date)) %>% 
  mutate(age_group = str_remove_all(age_group, "A")) %>% 
  mutate(age_group = case_when(
    age_group == "unbekannt" ~ NA_character_,
    age_group == "80+" ~ "80-99",
    TRUE ~ age_group
  )) %>% 
  mutate(gender = case_when(
    gender == "W" ~ "F",
    gender == "unbekannt" ~ NA_character_,
    TRUE ~ gender
  )) %>% 
  group_by(state, county, age_group, gender, date) %>% 
  summarise(cases = sum(cases),
            deaths = sum(deaths),
            recovered = sum(recovered)) %>% 
  ungroup() %>% 
  filter(cases >= 0 & deaths >= 0) %>%
  filter(date < today()) %>% 
  mutate(state = str_replace_all(state, "ü", "ue")) %>% 
  mutate(state = str_replace_all(state, "ä", "ae")) %>% 
  mutate(state = str_replace_all(state, "ö", "oe")) %>% 
  mutate(state = str_replace_all(state, "ß", "ss")) %>% 
  mutate(county = str_replace_all(county, "ü", "ue")) %>% 
  mutate(county = str_replace_all(county, "ä", "ae")) %>% 
  mutate(county = str_replace_all(county, "ö", "oe")) %>% 
  mutate(county = str_replace_all(county, "ß", "ss")) %>% 
  mutate(county = str_remove(county, "\\(.+\\)")) %>% 
  mutate(county = str_trim(county)) %>% dplyr::filter(!is.na(age_group))%>% dplyr::filter(!is.na(gender))

shape_county <- st_read("covid19-tracking-germany/de_county.shp", quiet = TRUE) %>% rename(county = GEN) %>% 
  dplyr::select(county, BEZ, geometry) %>% 
  mutate(county = as.character(county)) %>% 
  mutate(county = str_replace_all(county, "ü", "ue")) %>% 
  mutate(county = str_replace_all(county, "ä", "ae")) %>% 
  mutate(county = str_replace_all(county, "ö", "oe")) %>% 
  mutate(county = str_replace_all(county, "ß", "ss")) %>% 
  mutate(county = str_remove(county, "\\(.+\\)")) %>% 
  mutate(county = str_trim(county)) %>% 
  mutate(BEZ = case_when(
    BEZ == "Kreis" ~ "LK",
    BEZ == "Landkreis" ~ "LK",
    BEZ == "Stadtkreis" ~ "SK",
    BEZ == "Kreisfreie Stadt" ~ "SK"
  )) %>% 
  unite(county, BEZ, county, sep = " ", remove = TRUE) %>% group_by(county) %>% summarize()

# split(.$county) %>% 
# lapply(st_union) %>% 
# do.call(c, .) %>% st_sf()

foo <- covid_de %>% 
  mutate(county = case_when(
    county == "Region Hannover" ~ "LK Region Hannover",
    county == "SK Muelheim a.d.Ruhr" ~ "SK Muelheim an der Ruhr",
    county == "StadtRegion Aachen" ~ "LK Staedteregion Aachen",
    county == "SK Offenbach" ~ "SK Offenbach am Main",
    county == "LK Bitburg-Pruem" ~ "LK Eifelkreis Bitburg-Pruem",
    county == "SK Landau i.d.Pfalz" ~ "SK Landau in der Pfalz",
    county == "SK Ludwigshafen" ~ "SK Ludwigshafen am Rhein",
    county == "SK Neustadt a.d.Weinstrasse" ~ "SK Neustadt an der Weinstrasse",
    county == "SK Freiburg i.Breisgau" ~ "SK Freiburg im Breisgau",
    county == "LK Landsberg a.Lech" ~ "LK Landsberg am Lech",
    county == "LK Muehldorf a.Inn" ~ "LK Muehldorf a. Inn",
    county == "LK Pfaffenhofen a.d.Ilm" ~ "LK Pfaffenhofen a.d. Ilm",
    county == "SK Weiden i.d.OPf." ~ "SK Weiden i.d. OPf.",
    county == "LK Neumarkt i.d.OPf." ~ "LK Neumarkt i.d. OPf.",
    county == "LK Neustadt a.d.Waldnaab" ~ "LK Neustadt a.d. Waldnaab",
    county == "LK Wunsiedel i.Fichtelgebirge" ~ "LK Wunsiedel i. Fichtelgebirge",
    county == "LK Neustadt a.d.Aisch-Bad Windsheim" ~ "LK Neustadt a.d. Aisch-Bad Windsheim",
    county == "LK Dillingen a.d.Donau" ~ "LK Dillingen a.d. Donau",
    county == "LK Stadtverband Saarbruecken" ~ "LK Regionalverband Saarbruecken",
    county == "LK Saar-Pfalz-Kreis" ~ "LK Saarpfalz-Kreis",
    county == "LK Sankt Wendel" ~ "LK St. Wendel",
    county == "SK Brandenburg a.d.Havel" ~ "SK Brandenburg an der Havel",
    str_detect(county, "Berlin") ~ "SK Berlin",
    TRUE ~ county
  ))  

stateAndCounty <- foo %>% dplyr::select(state,
                                        county) %>% unique()

foo=foo%>% group_by(county, date) %>% 
  summarise(
    cases = sum(cases),
    deaths = sum(deaths)) %>% 
  ungroup() %>% 
  complete(county, date, fill = list(cases = 0, deaths = 0)) %>% 
  group_by(county) %>% 
  mutate(cumul_cases = cumsum(cases),
         cumul_deaths = cumsum(deaths)) %>% 
  ungroup() %>%  dplyr::filter(date==max(date))

shape_county <- shape_county %>% dplyr::left_join(foo, by="county") %>% dplyr::left_join(stateAndCounty, by="county")


population_state <- read.csv("covid19-tracking-germany/demographics_de.csv")

# https://www.citypopulation.de/en/germany/admin/
population_de_table <- read.csv("population_de_table.csv") %>% dplyr::filter(str_detect(Status, "County")) %>% 
  mutate(Name = str_replace_all(Name, "ü", "ue")) %>% 
  mutate(Name = str_replace_all(Name, "ä", "ae")) %>% 
  mutate(Name = str_replace_all(Name, "ö", "oe")) %>% 
  mutate(Name = str_replace_all(Name, "ß", "ss")) %>% 
  mutate(Name=gsub("\\s*\\([^\\)]+\\)","",Name)) %>% transmute(Name=Name, population=Population.Estimate.2018.12.31)

population_de_table$Name[grep("Heidekreis", population_de_table$Name)] <- "Heidekreis"
population_de_table$Name[grep("Bitburg-Pruem", population_de_table$Name)] <- "Eifelkreis Bitburg-Pruem"

library(stringr)
shape_county$stateShort=str_sub(shape_county$county,1,2)
shape_county$county=str_sub(shape_county$county,4,-1)

shape_county <- shape_county %>% dplyr::left_join(population_de_table, by=c("county"="Name"))

rate_de <- covid_de %>% group_by(age_group,gender) %>% summarise(totalCase=sum(cases,na.rm = TRUE)) 
pop_agg_state <- population_state %>% group_by(age_group, gender) %>% summarise(sumPop=sum(population)) 
rate_de$incidenceRate=rate_de$totalCase/pop_agg_state$sumPop

new_population_state <-  population_state%>% group_by(state,age_group, gender) %>% summarise(sumPop=sum(population)) 
new_population_state$e <- new_population_state$sumPop*rate_de$incidenceRate
new_population_state <- new_population_state %>% group_by(state) %>% summarise(E_state=sum(e))

pop_state <- shape_county %>% dplyr::group_by(state) %>% summarise(pop_state=sum(population)) %>% as.data.frame()

shape_county <- shape_county %>% dplyr::left_join(new_population_state, by="state") %>% dplyr::left_join(pop_state,by="state")
shape_county <- shape_county %>% mutate(E=E_state*population/pop_state) %>% dplyr::select(-cases, -deaths)
shape_county <- shape_county %>% mutate(popDensity=population/st_area(shape_county))

# the following was used to get "final_shp", they are in the same order
tempShp <- st_transform(shape_county, "+proj=longlat +datum=WGS84") 
tempShp <- as(shape_county,"Spatial")
tempShp <- spTransform(tempShp, "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# create new model data
new_model_data <- Final_shp@data
new_model_data$ID <- 1:nrow(new_model_data)
new_model_data <- new_model_data %>% dplyr::select(county, no2:ID)
new_model_data <- cbind(new_model_data, tempShp@data %>% dplyr::select(-county))

# set priors
prior <- list(
  theta1 = list(
    prior = "pc.prec",
    param = c(1, 0.5)),
  theta2 = list(
    prior = "gaussian",
    param = c(0, 1.8))
)
# leroux model
new_spatialModel_SIR <- inla(cumul_cases ~ no2+pm25+so2+Temperature+scale(Benzene)+scale(Aresenic)+scale(Cadmium)+scale(Nickel)+scale(popDensity)+f(ID, model = 'besagproper2', graph = W, hyper = prior), E = E, family = 'poisson', data = new_model_data
                              ,control.compute = list(config=TRUE,dic = TRUE, waic = TRUE)
                              ,control.predictor = list(compute = TRUE))
# leroux model result
round(100*(exp(new_spatialModel_SIR$summary.fixed[-1,c("mean","0.025quant","0.975quant")])-1),2)





