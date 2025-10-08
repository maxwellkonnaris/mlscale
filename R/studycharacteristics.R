
datasets <- list(
  "Vandeputte2021" = list(
    sequencingtype = "16S rRNA",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = NULL,
    ID = "Sample_ID",
    PMID = "34795283"
  ),
  "CvandeVelde2022" = list(
    sequencingtype = "16S rRNA",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = "Prefix",
    ID = "Sample",
    PMID = "37938658"
  ),
  "Vandeputte2017" = list(
    sequencingtype = "16S rRNA",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = "Status",
    ID = "Sample",
    PMID = "29143816"
  ),
  "Vandeputte2017 Disease Cohort" = list(
    sequencingtype = "16S rRNA",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = "Status",
    ID = "Sample",
    PMID = "29143816"
  ),
  "Vandeputte2017 Study Cohort" = list(
    sequencingtype = "16S rRNA",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = "Enterotype",
    ID = "Sample",
    PMID = "29143816"
  ),
  "Pereira2023" = list(
    sequencingtype = "16S rRNA",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = c("condition", "time_hours", "medium", "replicate"),
    ID = "Sample",
    PMID = "39572788"
  ),
  
  "Krawczyk2022" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Tick",
    sampletype = "Whole Organism",
    covariates = c("Life stage", "Location"),
    ID = "Sample_name",
    PMID = "35927748"
  ),
  "Liao2021" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = c("VanA", "Consistency", "Pool", "Timepoint"),
    ID = "SampleID",
    PMID = "33654104"
  ),
  
  "Stammler2016" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Mice",
    sampletype = "Fecal",
    covariates = c("Treatment", "Dilution", "Background"),
    ID = "SampleID",
    PMID = "27329048"
  ),
  
  "Dreier2022" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Cheese",
    sampletype = "Cheese",
    covariates = NULL,
    ID = "Sample_name",
    PMID = "35130830"
  ),
  
  "GALAXY" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = "cohort",
    ID = "ID",
    PMID = "39541968"
  ),
  
  "MetaCardis" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = c("PatientGroup", "SMOKE", "Gender", "CENTER", "AGE", "cohort", "pa_work_2cl"),
    ID = "Sample",
    PMID = c("35177860", "39541968", "34880489")
  ),

  "Marotz2021" = list(
    sequencingtype = "16S rRNA",
    loadtype = c("Flow Cytometry", "qPCR"),
    organismtype = "Human",
    sampletype = "Saliva",
    covariates = c("Treatment_code", "sex", "host_age", "mouthwash_regularly", "processing"),
    ID = "SampleID",
    PMID = "33594005"
  ),
  
  "Vieira_Silva2019" = list(
    sequencingtype = "16S rRNA",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = c("Diagnosis", "Age", "Gender", "BMI"),
    ID = "ID",
    PMID = "31209308"
  ),
  
  "Contijoch2019" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype = "qPCR",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = c("experimental_timepoint", "sample_mass"),
    ID = "Sample Name",
    PMID = "30666957"
  ),

  "Contijoch2019 16S" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "multi-species",
    sampletype = "Fecal",
    covariates = c( "Diagnosis", "experimental_timepoint", "sample_mass_mg",  "Organism"),
    ID = "Sample Name",
    PMID = "30666957"
  ),

  "Tunsakul2024" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = c("environment"),
    ID = "Sample",
    PMID = "38650647"
  ),
  
  "Alessandri2024" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = c("Fecal", "Vaginal"),
    covariates = c("treatmentfromlibraryname", "env_broad_scale"),
    ID = "Sample",
    PMID = "39364592"
  ),

  "Alessandri2024_fecal" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = c("Fecal"),
    covariates = c("treatmentfromlibraryname"),
    ID = "Sample",
    PMID = "39364592"
  ),

  "Alessandri2024_vaginal" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = c("Vaginal"),
    covariates = c("treatmentfromlibraryname"),
    ID = "Sample",
    PMID = "39364592"
  ),
  
  "Maghini2023" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype = "qPCR",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = c("Condition", "Donor"),
    ID = "ID",
    PMID = "37106038"
  ),
  
  "Garcia_Martinez2024" = list(
    sequencingtype = "16S rRNA",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = c("population"),
    ID = "Sample",
    PMID = "39462312"
  ),
  
  "Sternes2024" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Human",
    sampletype = "Gastrointestinal mucosa biopsy",
    covariates = c("Sex", "AGE", "BMI", "PPI", "NC_total_score"),
    ID = "Sample_name",
    PMID = "39469457"
  ),
  
  "Rao2021" = list( 
    sequencingtype = "16S rRNA",
    loadtype = c("Flow Cytometry", "qPCR", "MK_spike"),
    organismtype = c("Human", "Mouse", "Mock"),
    sampletype = "Fecal",
    covariates = c("sample_type"), # can probably mine more from combinedphasesamplename
    ID = c("Sample_name"),
    PMID = "33627867"
  ),
  
  "Tettamanti_Boshier2020" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Human",
    sampletype = "Vaginal muscosa",
    covariates = c("Hours_In_Study"),
    ID = "Sample_ID",
    PMID = "32265316"
  ),
  
  "Kruger2024" = list(
    sequencingtype = "16S rRNA",
    loadtype = "ddPCR",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = c("BristolStoolScale_highest","Timepoint", "pH"),
    ID = "SampleID",
    PMID = "39427011"
  ),
  
  "Liu2017" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Human",
    sampletype = "Coronal sculcus",
    covariates = c("HIV Status"),
    ID = "Sample_name",
    PMID = "28743816"
  ),
  
  "Fu2023" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype = "Flow Cytometry",
    organismtype = "Mock", 
    sampletype = "Wastewater",
    covariates = NULL,
    ID = "Sample",
    PMID = "38868349"
  ),
  
  "Jin2022" = list(
    sequencingtype = "16S rRNA",
    loadtype = "ddPCR",
    organismtype = "Mice",
    sampletype = "Cecum content",
    covariates = c("diet","location", "sample_type", "subject"),
    ID = "Sample",
    PMID = "35194029"
  ),
  
  "Zaramela2022" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = "Saliva",
    covariates = c("gender", "Pool"),
    ID = "ID",
    PMID = "36317886"
  ),
  
  "Feng2023" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Chicken",
    sampletype = "Gut segment",
    covariates = c("Segment", "Date"),
    ID = "Sample",
    PMID = "38868437"
  ),
  
  "Reese2022" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Chimpanzee",
    sampletype = "Fecal",
    covariates = c("Age group", "Community", "Sex", "YearQuarter"),
    ID = "Sample_name",
    PMID = "33232664"
  ),
  
  "Barlow2020" = list(
    sequencingtype = "16S rRNA",
    loadtype = "ddPCR",
    organismtype = "Mice",
    sampletype = "Fecal",
    covariates = c("Diet", "Site"),
    ID = "Sample",
    PMID = "32444602"
  ),

  "Barlow2021" = list(
    sequencingtype = "16S rRNA",
    loadtype = "ddPCR",
    organismtype = "Human",
    sampletype = c("Duodenal-Aspirate","Duodenum-Saliva"),
    covariates = NULL,
    ID = "Sample",
    PMID = "34724979"
  ),

  "Morton2019" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype = c("qPCR", "Flow Cytometry"),
    organismtype = "Human",
    sampletype = "Saliva",
    covariates = c("processing", "sex", "Host_age", "mouthwash_regularly"),
    ID = "Accession",
    PMID = "31222023"
  ),

  "Morton2019 16S" = list(
    sequencingtype = "16S rRNA",
    loadtype = c("qPCR", "Flow Cytometry"),
    organismtype = "Human",
    sampletype = "Saliva",
    covariates = c("processing", "sex", "Host_age", "mouthwash_regularly"),
    ID = c("Accession"),
    PMID = "31222023"
  ),
  
  "Prochazkova2024" = list(
    sequencingtype = "16S rRNA",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = c("Fecal"),
    covariates = c("Bristol_scale_Mean", "Faecal_pH_Mean", "Stool_freq_Mean","Stool_moisture_Mean"), # We dont have access to any of the other metadata..
    ID = "ID",
    PMID = "39604623"
  ),
  
  "Zemb2020" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Pig",
    sampletype = "Fecal",
    covariates = c("mgfeces","added_Coli"), # We dont have access to any of the other metadata..
    ID = "Sample.Name",
    PMID = "31927795"
  ),
  
  "Jin2024" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Human",
    sampletype = "Semen",
    covariates = c("Fertility", "Diagnosis", "Group", "Sperm_concentration", "Sperm_motility", "Age"),
    ID = "Sample",
    PMID = "39359400"
  ),
  
  "Galazzo2020" = list(
    sequencingtype = "16S rRNA",
    loadtype = c("qPCR", "ddPCR", "Flow Cytometry"),
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = c("sample_treatment", "host_subject_id"),
    ID = "Sample_name",
    PMID = "32850498"
  ),
  
  "Lin2019" = list(
    sequencingtype = "16S rRNA",
    loadtype = "Flow Cytometry",
    organismtype = "Marine",
    sampletype = "Seawater",
    covariates = c("Station", "Line","Filtered Seawater Vol [L]", "SurfChl [mg m-3]"),
    ID = "SampleID",
    PMID = "30552195"
  ),
  
  "Suriano2022" = list(
    sequencingtype = "16S rRNA",
    loadtype = "Flow Cytometry",
    organismtype = "Mice",
    sampletype = "Fecal",
    covariates = c("Day"),
    ID = "Sample_name",
    PMID = "36516223"
  ),
  
  "Thiruppathy2025" = list(
    sequencingtype = "Shotgun Metagenomics",
    loadtype = "Flow Cytometry",
    organismtype = "Human",
    sampletype = "Skin",
    covariates = c("PMA_treated", "BodySite_ID", "Age", "Sex", "Exercise_per_week", "Face_sunscreen", "Body_shower_fq_by_day"),
    ID = "Sample",
    PMID = "40038838"
  ),
  
  "Wagner2025" = list(
    sequencingtype = "16S rRNA",
    loadtype = "Flow Cytometry",
    organismtype = "Piglet",
    sampletype = "Fecal",
    covariates = NULL,
    ID = "Sample",
    PMID = "40196033"
  ),

  "Rolling2021" = list(
    sequencingtype = "16S rRNA",
    loadtype = "qPCR",
    organismtype = "Human",
    sampletype = "Fecal",
    covariates = NULL,
    ID = "Sample",
    PMID = "34764444"
  ),
  
  "Kallastu2023" = list(
    sequencingtype = "16S rRNA",
    loadtype = c("Flow Cytometry", "qPCR", "CFU"),
    organismtype = "Food",
    sampletype = "Food",
    covariates = c("pma_treatment", "Organism"),
    ID = "Sample",
    PMID = "36691592"
  )
)
