

# Split dataset into  subtrials
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

load("datalist.RData") #specify path if needed
datalist_full <- datalist
for (j in 1:datalist_full$nsub) {
  x = datalist_full$x[datalist_full$sub == j]
  y = datalist_full$y[datalist_full$sub == j]
  b = datalist_full$b[datalist_full$sub == j]
  t = datalist_full$t[datalist_full$sub == j]
  N = datalist_full$N/datalist_full$nsub
  resp = datalist_full$resp[datalist_full$sub == j]
  sub = datalist_full$sub[datalist_full$sub == j]
  nsub = 1
  rr_c = datalist_full$rr_c[j]
  rr_t = datalist_full$rr_t[j]
  datalist <- list(x = x, y = y, b = b, t = t, N = N, resp = resp, sub = sub, nsub = nsub, rr_c = rr_c, rr_t = rr_t)
  save(datalist, file = paste("datalist", j, ".RData", sep = ""))
}
