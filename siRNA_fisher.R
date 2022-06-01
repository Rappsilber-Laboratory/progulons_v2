
test_matrix <- matrix( c(12, 10, 3, 21), 
                       nrow = 2, 
                       dimnames = list( c("siRNA_high", "siRNA_low"),
                                        c("prn_high", "prn_low")))
test_matrix
fisher.test(x = test_matrix, alternative = "greater")
