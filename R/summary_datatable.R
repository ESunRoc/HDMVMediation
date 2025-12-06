to_table <- function(mod_boot_summ, mod_num, p){
  mod_boot_summ <- as.data.frame(mod_boot_summ)
  
  bca_inter <- paste0("(",round(mod_boot_summ$bca_lowerCL,4),", ",round(mod_boot_summ$bca_upperCL,4),")")
  piv_inter <- paste0("(",round(mod_boot_summ$piv_lowerCL,4),", ",round(mod_boot_summ$piv_upperCL,4),")")
  
  if(mod_num %in% c(1,2)) {
    mod_boot_outcome <- c(rep("PI", p+2), rep("ICDAS", p+2), rep("DS", p+2), rep("MS", p+2), rep("FS", p+2))
  } else if(mod_num %in% c(3,4,5)) {
    mod_boot_outcome <- c(rep("ABO", p+2), rep("BW", p+2))
  } else{
    mod_boot_outcome <- c(rep("PI", p+2), rep("ICDAS", p+2), rep("DS", p+2), rep("MS", p+2), rep("FS", p+2), rep("ABO", p+2), rep("BW", p+2))
  }
  
  mod_boot_table_df <- data.frame("Outcome" = mod_boot_outcome,
                                  "Estimand" = gsub("(_resp[1-9])", "", gsub("(_ide_resp[1-9])", "", rownames(mod_boot_summ))),
                                  "Orig. Est." = round(mod_boot_summ$OrigEst,4),
                                  "Mean(boot)" = round(mod_boot_summ$Mean_boot,4),
                                  "sd(boot)" = round(mod_boot_summ$boot_SE,4),
                                  "pval" = signif(mod_boot_summ$boot_pval,4),
                                  "BC_a CI" = bca_inter,
                                  "PBCI" = piv_inter)
  mod_boot_table_df_print <- mod_boot_table_df[which(mod_boot_table_df$Orig..Est.!=0),]
  
  
  DT::datatable(mod_boot_table_df_print, rownames = F,
                colnames = c("Outcome", "Estimand", "Orig. Est.", 
                             "Mean(boot)", "sd(boot)", "p-value",
                             "BCa CI", "PBCI"))
  # clipr::write_clip(knitr::kable(mod_boot_table_df_print, digits = 4, format = "latex", align = "c"))
}