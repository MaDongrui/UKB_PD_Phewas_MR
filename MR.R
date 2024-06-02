library(TwoSampleMR)
library(MRPRESSO)
library(ggplot2)

# Define the list of datasets to process
datasets <- c(
  'ukb-b-5238', 'ukb-b-1431', 'ukb-b-4650', 'ukb-b-6306', 'ukb-b-18335', 'ukb-b-10387',
  'ukb-b-2091', 'ukb-b-2750', 'ukb-b-918', 'ukb-b-13352', 'ukb-b-17241', 'ukb-b-11085',
  'ukb-b-298', 'ukb-b-2134', 'ukb-b-6244', 'ukb-b-15926', 'ukb-b-8121', 'ukb-b-5779',
  'ukb-b-5239', 'ukb-b-5716', 'ukb-b-16878', 'ukb-b-3460', 'ukb-b-4256', 'ukb-b-2047',
  'ukb-b-4801', 'ukb-b-3831', 'ukb-b-14180', 'ukb-b-13745', 'ukb-b-19809', 'ukb-b-8476',
  'ukb-b-10169', 'ukb-b-1419', 'ukb-b-18674', 'ukb-b-8196', 'ubm-a-2715', 'ubm-a-2716',
  'ubm-a-2717', 'ubm-a-2718', 'ubm-a-2719', 'ubm-a-2720', 'ubm-a-2721', 'ubm-a-2722',
  'ubm-a-2723', 'ubm-a-2724', 'ubm-a-2725', 'ubm-a-2726', 'ubm-a-2727', 'ubm-a-2728',
  'ubm-a-2729', 'ubm-a-2730', 'ubm-a-2731', 'ubm-a-2732', 'ubm-a-2734', 'ubm-a-2735',
  'ubm-a-2736', 'ubm-a-2737', 'ubm-a-2738', 'ubm-a-2739', 'ubm-a-2740', 'ubm-a-2741',
  'ubm-a-2742', 'ubm-a-2743', 'ubm-a-2744', 'ubm-a-2745', 'ubm-a-2931', 'ubm-a-2932',
  'ubm-a-2941', 'ubm-a-2950', 'ubm-a-2951', 'ubm-a-2952', 'ubm-a-2955', 'ubm-a-2821',
  'ubm-a-2822', 'ubm-a-2823', 'ubm-a-2825', 'ubm-a-2826', 'ubm-a-2827', 'ubm-a-2828',
  'ubm-a-2829', 'ubm-a-2830', 'ubm-a-2831', 'ubm-a-2832', 'ubm-a-2833', 'ubm-a-2834',
  'ubm-a-2835', 'ubm-a-2836', 'ubm-a-2837', 'ubm-a-2838', 'ubm-a-2839', 'ubm-a-2840',
  'ubm-a-2841', 'ubm-a-2842', 'ubm-a-2843', 'ubm-a-2844', 'ubm-a-2845', 'ubm-a-2846',
  'ubm-a-2847', 'ubm-a-2848', 'ubm-a-2849', 'ubm-a-2850', 'ubm-a-2851', 'ubm-a-3037',
  'ubm-a-3056', 'ukb-b-10787', 'ukb-b-20124', 'ukb-b-7953', 'ukb-b-4080', 'ukb-b-19234',
  'ukb-b-8875', 'ukb-b-7992', 'ukb-b-2122', 'ukb-b-8607', 'ukb-b-13354', 'ukb-b-14540',
  'ukb-b-7376', 'ukb-b-14068', 'ukb-b-12828', 'ukb-b-14310', 'ukb-b-16099', 'ukb-b-17271',
  'ukb-b-19520', 'ukb-b-16698', 'ukb-b-19925', 'ukb-b-9093', 'ukb-b-17409', 'ukb-b-9685',
  'ukb-d-30000_irnt', 'ukb-d-30020_irnt', 'ukb-d-30030_irnt', 'ukb-d-30140_irnt',
  'ukb-d-30200_irnt', 'ukb-d-30240_irnt', 'ukb-d-30250_irnt', 'ukb-d-30280_irnt',
  'ukb-d-30290_irnt', 'ukb-a-333', 'ukb-a-335', 'ukb-b-8951', 'ukb-d-30810_irnt',
  'ukb-d-30880_irnt', 'ukb-b-5945', 'ukb-b-419', 'ukb-b-185', 'ukb-b-6134', 'ukb-b-188',
  'ukb-b-16489', 'ukb-b-7926', 'ukb-b-17729', 'ukb-b-13824', 'ukb-b-16671'
)

# Create a file to record error messages
error_log <- file("error_log.txt", open = "w")

# Process each dataset
for (dataset in datasets) {
  cat(paste("Processing dataset: ", dataset, "\n"))
  
  tryCatch({
    # Processing code, replacing the dataset name with `dataset`
    
    ukb_b <- extract_instruments(
      outcomes = dataset,  # Use the current dataset name
      clump = TRUE, r2 = 0.001, p1 = 1e-6,
      kb = 10000, access_token = NULL
    )
    
    ukb_b$R2 <- 2 * (1 - ukb_b$eaf.exposure) * ukb_b$eaf.exposure * (ukb_b$beta.exposure)^2
    ukb_b$F <- (ukb_b$R2) / (1 - ukb_b$R2) * (ukb_b$samplesize.exposure - 2)
    
    # Ensure F > 10
    ukb_b <- subset(ukb_b, F > 10)
    
    write.csv(ukb_b, file = paste0(dataset, "_R2_F.csv"))
    
    pd_out <- extract_outcome_data(
      snps = ukb_b$SNP,
      outcomes = 'ieu-b-7',
      proxies = FALSE,
      maf_threshold = 0.01,
      access_token = NULL
    )
    
    mydata <- harmonise_data(
      exposure_dat = ukb_b,
      outcome_dat = pd_out,
      action = 2
    )
    
    res <- mr(mydata)
    generate_odds_ratios(res)
    mr_res <- generate_odds_ratios(res)
    write.csv(mr_res, file = paste0(dataset, "→pd.csv"))
    
    # Heterogeneity test
    het <- mr_heterogeneity(mydata)
    print(het)
    write.csv(het, file = paste0(dataset, "_heterogeneity.csv"))
    pleio <- mr_pleiotropy_test(mydata)
    pleio
    write.csv(pleio, file = paste0(dataset, "_pleiotropy.csv"))
    
    # Calculate SNP-exposure correlation
    mydata$r.exposure <- get_r_from_bsen(mydata$beta.exposure, mydata$se.exposure, mydata$samplesize.exposure)
    
    # Calculate SNP-outcome correlation
    mydata$r.outcome <- get_r_from_bsen(mydata$beta.outcome, mydata$se.outcome, mydata$samplesize.outcome)
    
    # Run directionality_test using the calculated correlations
    Steiger <- directionality_test(mydata)
    write.csv(Steiger, file = paste0(dataset, "_Steiger.csv"))
    
    # Perform mr_presso function
    mr_presso <- capture.output({
      mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure",
                OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = mydata, NbDistribution = 1000,
                SignifThreshold = 0.05)
    })
    
    # Specify the file path to save the output
    output_file <- paste0(dataset, "_mr_presso.txt")
    
    # Write the output to a file
    writeLines(mr_presso, con = output_file)
    
    p1 <- mr_scatter_plot(res, mydata)
    ggsave(p1[[1]], file = paste0(dataset, "_scatter_plot.pdf"), width = 7, height = 7)
 
}, error = function(e) {
    # If an error occurs, record the error message in the log file
    cat(paste("Error processing dataset: ", dataset, "\n"))
    cat(paste("Error message: ", conditionMessage(e), "\n"))
    cat("\n", file = error_log, append = TRUE)
  })
}

# Close the error log file
close(error_log)
