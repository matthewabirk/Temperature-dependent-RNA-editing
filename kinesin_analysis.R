khc_vel = readxl::read_xlsx('Data/kinesin_data.xlsx', sheet = 'Velocity', skip = 1)

t.test(x = khc_vel$`WT-RT`, y = khc_vel$`K282R-RT`, var.equal = TRUE)
t.test(x = khc_vel$`WT-11C`, y = khc_vel$`K282R-11C`, var.equal = TRUE)

t.test(x = khc_vel$`WT-RT`, y = khc_vel$`WT-11C`, var.equal = TRUE)
t.test(x = khc_vel$`K282R-RT`, y = khc_vel$`K282R-11C`, var.equal = TRUE)



khc_rl = readxl::read_xlsx('Data/kinesin_data.xlsx', sheet = 'Run length ', skip = 1)

t.test(x = khc_rl$`WT-RT`, y = khc_rl$`K282R-RT`, var.equal = TRUE)
t.test(x = khc_rl$`WT-11C`, y = khc_rl$`K282R-11C`, var.equal = TRUE)

t.test(x = khc_rl$`WT-RT`, y = khc_rl$`WT-11C`, var.equal = TRUE)
t.test(x = khc_rl$`K282R-RT`, y = khc_rl$`K282R-11C`, var.equal = TRUE)
