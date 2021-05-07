library(sfcr)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(kableExtra)
library(grid)

grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}




eqs <- sfcr_set(
  #model creation,
  #(1) arterial flows,
  #income,
  #yus~ cus+ gus +xus- imus,
  #yla ~ cla+ gla+ xla -imla,
  #disposable income,
  ydus~(yus + fbus +rus[-1]*busus_d[-1]+ rdus[-1]*dus_d[-1] +cba_d[-1]+cgus)*(1-thetaus),
  ydla~(yla + fbla + rla[-1]*blala_d[-1] + rdla[-1]*dla_d[-1]+der_d[-1]+cgla)*(1-thetala), 
  #disposable income (hs),
  #ydus_hs~ydus+ d(xrla)*busla_s[-1],
  #ydla_hs~ydla+ d(xrus)*blaus_s[-1],
  #wealth,
  vus ~ vus[-1]+  ydus-cus,
  vla ~ vla[-1]+ ydla-cla,
  #capital gains,
  cgus ~ (pcba-pcba[-1])*cba_s[-1],
  cgla ~ (pder-pder[-1])*der_s[-1],
  #tax,
  tla ~ thetala*(yla+ fbla+rla[-1]*blala_d[-1]+ rdla[-1]*dla_d[-1]+der_d[-1]),
  tus ~ thetaus*(yus+fbus +rus[-1]*busus_d[-1]+rdus[-1]*dus_d[-1] + cba_d[-1]) ,
  #cb profits,
  fcbla ~ rla[-1]*bcbla_d[-1] + rus[-1]*bcblaus_s[-1]*xrus,
  fcbus ~ rus[-1]*bcbus_d[-1], 
  #government budget constraint,
  bla_s ~ bla_s[-1]+gla-tla+ rla[-1]*bla_s[-1] - fcbla, 
  bus_s ~ bus_s[-1]+gus-tus+ rus[-1]*bus_s[-1] - fcbus,
  #current and capital account,
  cabla~xla - imla - rla[-1]*bbusla_s[-1] + rus[-1]*bblaus_s[-1]*xrus + rus[-1]*bcblaus_s[-1]*xrus ,
  kabla ~ (bbusla_s-bbusla_s[-1]) - (bcblaus_s-bcblaus_s[-1])*xrus- (bblaus_s-bblaus_s[-1])*xrus,
  cabus~xus - imus + rla[-1]*bbusla_s[-1]*xrla - rus[-1]*bblaus_s[-1] - rus[-1]*bcblaus_s[-1],
  #kabus ~ -d(bbla_s)*xrla + d(blaus_s) + d(bcblaus_s),
  #(2) trade,
  #export and import,
  pmla~ exp(mi_0la+  mi_1la *log(r_pusy)  + (1- mi_1la)*log(r_play) - mi_1la *log(xrla)),
  pxla ~exp(chi_0la+  chi_1la *log(r_pusy)  + (1- chi_1la)*log(r_play)- chi_1la *log(xrla)),
  pxus ~ pmla*xrla,  
  pmus ~ pxla*xrla,
  r_xla~ exp(e_la - eta_la *log(pmus/r_pusy) + epsilon_la*log(r_yus)),
  r_imla ~ exp(p_la - psi_la *log(pmla[-1]/r_play[-1]) + pi_la*log(r_yla)),
  r_xus  ~r_imla,
  r_imus~r_xla,
  xla ~ r_xla*pxla,
  xus ~ r_xus*pxus,
  imla ~ r_imla*pmla,  
  imus ~ r_imus*pmus,
  #(3) income and expenditure,
  #real disposable income is of the haig-simon type,
  r_vla ~ vla/r_pdsla, 
  r_vus~ vus/r_pdsus,
  r_ydla ~ ydla/r_pdsla - r_vla[-1]* (r_pdsla-r_pdsla[-1])/r_pdsla ,
  r_ydus ~ ydus/r_pdsus - r_vus[-1]* (r_pdsus-r_pdsus[-1])/r_pdsus ,
  r_cla ~ alpha_1la*r_ydlae + alpha_2la*r_vla[-1],
  r_cus ~alpha_1us*r_yduse + alpha_2us*r_vus[-1],
  r_ydlae~ (r_ydla + r_ydla[-1])/2,
  r_yduse~ (r_ydus + r_ydus[-1])/2,
  r_sla ~r_cla+r_gla+r_xla ,
  r_sus ~r_cus+r_gus+r_xus ,
  sla ~r_sla*r_plas,
  sus ~r_sus*r_puss,
  r_plas ~((1+phila)*(wla*nla+imla)) /r_sla,
  r_puss ~((1+phius)*(wus*nus+imus)) / r_sus,
  r_pdsla ~(sla-xla) / (r_sla-r_xla),
  r_pdsus ~(sus-xus) / (r_sus-r_xus),
  dsla ~sla-xla ,
  dsus ~sus-xus ,
  r_dsla ~r_cla+r_gla ,
  r_dsus ~r_cus+r_gus  ,
  yla~sla-imla,
  yus ~sus-imus,
  r_yla ~r_sla-r_imla ,
  r_yus ~r_sus-r_imus ,
  r_play ~ yla /r_yla ,
  r_pusy ~yus/ r_yus,
  cla ~r_cla*r_pdsla ,
  cus ~r_cus*r_pdsus  ,
  gla ~r_gla*r_pdsla ,
  gus ~r_gus*r_pdsus ,
  nla ~ r_yla/ r_prla , 
  nus ~ r_yus/ r_prus ,
  #(4) financial intermediaries,
  #us financial sector,
  bbus_d~dus_s-bbusla_d-hbus_d+pcba*cba_s,
  #d(vbus)~fbus,
  bbusla_d~rho_2us*cba_s,
  hbus_d~rho_0us*dus_s,
  rdus~ rdus[-1] + rho_1us*(rus-rus[-1]),
  cgbus ~ (xrla-xrla[-1])*bbusla_s[-1],
  pcba~(1/rcba)+log(pmus),
  #d(rcba)~rho_0int*,
  #d(rlus)~rho_2us*d(rus),
  fbus~rus[-1]*bbus_d[-1]-rdus*dus_s[-1]+rla[-1]*bbusla_d[-1]*xrla +cgbus-cba_d[-1]-cgus,
  #la financial sector,
  bbla_d~dla_s-bblaus_d-hbla_d+pder*der_s,
  #d(vbla)~fbla,
  bblaus_d~rho_2la*der_s,
  hbla_d~rho_0la*dla_s,
  rdla ~ rdla[-1]+ rho_1la*(rla-rla[-1]),
  cgbla ~ (xrus-xrus[-1])*bblaus_s[-1],
  pder~1/rder,
  #d(rder)~-rho_2int*(rla-rus),
  #d(rlla)~rho_2la*d(rla),
  fbla~rla[-1]*bbla_d[-1]-rdla*dla_s[-1]+rus[-1]*bblaus_d[-1]*xrus +cgbla-der_d[-1]-cgla,
  #intermediaries,
  #d(vint)~rus[-1]*bintus_d[-1]+rla[-1]*bintla_d[-1]*xrla +cgint-rcba[-1]*cba_d[-1],
  #cba_s~cba_d,
  #cgint~d(xrla)*bintla_s[-1],
  #pcba~1/rcba,
  #d(rcba)~rho_0int*d(pxla)-rho_1int*(rla-rus),
  #fint~rus[-1]*bintus_d[-1]+rla[-1]*bintla_d[-1]*xrla +cgint-rcba[-1]*cba_d[-1],
  #bintla_d~cba_s-bintus_d,
  #(5) assets demand,
  #asset demand for la resident,
  blala_d~vla*(lambda_10la+lambda_11la*rla-lambda_13la*rdla-lambda_14la*(rder+dxrlae)),
  dla_d~vla*(lambda_40la-lambda_41la*rla+lambda_43la*rdla-lambda_44la*(rder+dxrlae)),
  der_d~(vla/pder)*(lambda_50la-lambda_51la*rla-lambda_53la*rdla+lambda_54la*(rder+dxrlae)),
  #blaus_d~vla*(lambda_20la-lambda_21la*rla-lambda_23la*rdla+lambda_22la*(rus+dxruse)),
  #hla_d/ vla~lambda_30la - lambda_31la*(rus+dxruse) -lambda_32la*rla,
  #asset demand for us resident,
  busus_d~vus*(lambda_10us +lambda_11us*rus-lambda_13us*rdus-lambda_14us*(rcba+dxruse)),
  dus_d~vus*(lambda_40us-lambda_41us*rus+lambda_43us*rdus-lambda_44us*(rcba+dxruse)),
  cba_d~(vus/pcba)*(lambda_50us-lambda_51us*rus-lambda_53us*rdus+lambda_54us*(rcba+dxruse)),
  #busla_d~vus*(lambda_20us -lambda_21us*rus+lambda_22us*(rla+dxrlae)),
  #hus_d/ vus~lambda_30us -lambda_31us*rus -lambda_32us* (rla+dxrlae),
  #asset demand for intermediaries,
  #bintus_d~vint*(lambda_10int + lambda_11int*rus - lambda_12int*(rla+dxrlae)),
  #bintla_d~vint*(lambda_20int - lambda_21int*rus + lambda_22int*(rla+dxrlae)),
  #expected change in the exchange rate (expectations to update),
  dxrlae~d(pder)/ pder,
  dxruse~d(pcba)/ pcba,
  #(6) assets supply,
  #demand for cash,
  hus_d~ vus - busus_d- dus_d- pcba*cba_d ,
  hla_d~ vla - blala_d-dla_d -pder*der_d,
  #cb demand for b, h; d and cba,
  hus_s~ hus_d,
  hla_s~ hla_d,
  dus_s~dus_d,
  dla_s~dla_d,
  busus_s~ busus_d,
  blala_s ~blala_d,
  cba_s~cba_d,
  der_s~der_d,
  bbus_s~ bbus_d,
  bbla_s ~bbla_d,
  #aus_s ~aus_d,
  #ala_s ~ala_d,
  #bintus_s~bintus_d,
  #lus_s~lus_d,
  #lla_s~lla_d,
  hbus_s~hbus_d,
  hbla_s~hbla_d,
  #supply of domestic t bills to cb,
  bcbus_s~ bcbus_d,
  bcbla_s~ bcbla_d ,
  bcbus_d ~ bcbus_d[-1]+ (hus_s-hus_s[-1])+(hbus_s-hbus_s[-1]),
  bcbla_d ~ bcbla_d[-1]+ (hla_s-hla_s[-1])+(hbla_s-hbla_s[-1])-(bcblaus_s-bcblaus_s[-1])*xrus,
  #supply of assets abroad ,
  #busla_s~ bla_s- blala_s- bcbla_s,
  #exchange rate,
  #xrus~ blaus_d /blaus_s,
  xrla~ 1/xrus,
  bblaus_s~bblaus_d*xrla,
  bcblaus_d~bcblaus_s*xrus,
  xrus~ (bbusla_s)/(bbusla_d),
  bbusla_s~ bla_s - blala_s- bcbla_s -  bbla_s,
  #final equation (not possible to appear twice),
  #blaus_s~ blaus_d /xrus,
  rerla~(r_play/r_pusy)*(xrla),
  rerus~1/rerla,
  tbla~xla - imla,
  tbus~xus - imus,
  psbrla~d(bla_s),
  psbrus~d(bus_s),
  prbla~cabla+d(bla_s),
  prbus~cabus+d(bus_s),
  nwla~((vla)-(bla_s) +(bbla_d-dla_s+bblaus_d+hbla_d-pder*der_s)+(bcblaus_d+bcbla_d-hla_s)),
  nwus~((vus)-(bus_s) +(bbus_d-dus_s+bbusla_d+hbus_d-pcba*cba_s)+(bcbus_d-hus_s)),
  rsh_cla~r_cla/yla,
  sh_cabla~cabla/yla,
  sh_tbla~tbla/yla,
  sh_govdef~-psbrla/yla,
  sh_prbla~prbla/yla,
  nip~rla[-1]*bla_s[-1]/yla,
  govdeb~bla_s/yla,
  sh_bbusla_s~bbusla_s/bla_s,
  sh_cba~pcba*cba_d/r_vus,
  totla~pxla/pmla
)

external <- sfcr_set(
  alpha_1la ~ 0.9,
  alpha_1us ~ 0.6,
  alpha_2la ~ 0.053333276,
  alpha_2us ~ 0.213326723889398,
  #vertical and horizontal constraints respected
  lambda_10la ~ 0.4,
  lambda_11la ~ 0.5,
  lambda_12la ~ 0.5,
  lambda_13la~ 0.25,
  lambda_14la~ 0.25,
  lambda_20la ~ 0.3,
  lambda_21la ~ 0.5,
  lambda_22la ~ 0.5,
  lambda_23la~ 0.5,
  lambda_40la ~ 0.4,
  lambda_41la ~ 0.25,
  lambda_42la ~ 5,
  lambda_43la~ 0.5,
  lambda_44la~ 0.25,
  lambda_50la ~ 0.1,
  lambda_51la ~ 0.25,
  lambda_52la ~ 5,
  lambda_53la~ 0.25,
  lambda_54la~ 0.5,
  #vertical and horizontal constraints respected
  lambda_24la~ 5,
  lambda_30la ~ 0.25,
  lambda_31la ~ 5,
  lambda_32la ~ 5,
  lambda_33la~ 5,
  lambda_34la~ 5,
  #vertical and horizontal constraints respected
  lambda_10us ~ 0.4,
  lambda_11us ~ 0.5,
  lambda_12us ~ 0.5,
  lambda_13us~ 0.25,
  lambda_14us~ 0.25,
  lambda_20us ~ 0.3,
  lambda_21us ~ 0.5,
  lambda_22us ~ 0.5,
  lambda_23us~ 0.5,
  lambda_40us ~ 0.4,
  lambda_41us ~ 0.25,
  lambda_42us ~ 5,
  lambda_43us~ 0.5,
  lambda_44us~ 0.25,
  lambda_50us ~ 0.1,
  lambda_51us ~ 0.25,
  lambda_52us ~ 5,
  lambda_53us~ 0.25,
  lambda_54us~ 0.5,
  #vertical and horizontal constraints respected
  lambda_24us~ 5,
  lambda_30us ~ 0.25,
  lambda_31us ~ 5,
  lambda_32us ~ 5,
  lambda_33us~ 5,
  lambda_34us~ 5,
  #vertical and horizontal constraints respected
  lambda_60us~0.5,
  lambda_61us~0.5,
  lambda_62us~0.5,
  lambda_70us~0.5,
  lambda_71us~0.5,
  lambda_72us~0.5,
  lambda_60la~0.5,
  lambda_61la~0.5,
  lambda_62la~0.5,
  lambda_70la~0.5,
  lambda_71la~0.5,
  lambda_72la~0.5,
  #vertical and horizontal constraints respected
  e_la ~ - 1.18,
  eta_la ~ 0.7,
  epsilon_la ~ 0.8,
  p_la ~ - 3.015,
  psi_la ~ 0.7,
  pi_la ~  1.2,
  mi_0la ~ - 0.00001,
  chi_0la ~ - 0.00001,
  mi_1la ~ 0.7,
  chi_1la ~ 0.2,
  phila ~ 0.2381,
  phius ~ 0.2381,
  thetala ~ 0.2,
  thetaus ~ 0.2,
  thetabla ~ 0,
  thetabus ~ 0,
  rho_0us~0.03,
  rho_1us~0.9,
  rho_2us~2.7,
  rho_0la~0.03,
  rho_1la~0.9,
  rho_2la~2.7,
  rho_0int~0.9,
  rho_1int~0.9,
  rho_2int~0.9,
  # exogenous variables
  bcblaus_s ~ 0.02031,
  r_gla ~ 16,
  r_gus ~ 16,
  r_prla ~ 1.3333,
  r_prus ~ 1.3333,
  rla ~ 0.03,
  rus ~ 0.03,
  rcba~0.16,
  rder~0.16,
  wla ~ 1,
  wus ~ 1
)



initial <- sfcr_set(
  # starting values for stocks
  bcblaus_d ~ 0.02031,
  bla_s ~ 145.954295,
  blala_d ~ 102.18,
  blala_s ~ 102.18,
  bbus_d ~ 13.250255,
  bbus_s ~ 13.250255,
  bus_s ~ 146.005715,
  busus_d ~ 102.19,
  busus_s ~ 102.19,
  bcbla_d ~ 17.27839,
  bcbla_s ~ 17.27839,
  bcbus_d ~ 17.2995,
  bcbus_s ~ 17.2995,
  hla_d ~ 7.2987,
  hla_s ~ 7.2987,
  hus_d ~ 7.2995,
  hus_s ~ 7.2995,
  bblaus_d ~ 13.24565,
  bblaus_s ~ 13.24565,
  bbusla_d ~ 13.250255,
  bbusla_s ~ 13.250255,
  cba_d ~ 3.9178,
  cba_s ~ 3.9178,
  der_d ~ 3.9178,
  der_s ~ 3.9178,
  bbla_d ~ 13.24565,
  bbla_s ~ 13.24565,
  hbus_d~10,
  hbus_s~10,
  hbla_d~10,
  hbla_s~10,
  dus_s~12.01426,
  dla_s~12.01426,
  dus_d~12.01426,
  dla_d~12.01426,
  r_vla ~ 152.62,
  r_vus ~ 152.63,
  vla ~ 145.97921,
  vus ~ 145.99001,
  # other endogenous
  r_cla ~ 81.393,
  r_cus ~ 81.401,
  cabla ~ 0,
  cabus ~ 0,
  cla ~ 77.851,
  cus ~ 77.86,
  r_dsla ~ 97.393,
  r_dsus ~ 97.401,
  dsla ~ 93.154,
  dsus ~ 93.164,
  #dxrlae ~ 0,
  fcbla ~ 0.00869,
  fcbus ~ 0.00895,
  fbus~0.6,
  fbla~0.8,
  gla ~ 15.304,
  gus ~ 15.304,
  r_imla ~ 11.928,
  r_imus ~ 11.926,
  imla ~ 11.407,
  imus ~ 11.409,
  kabla ~ 0.00002,
  #kabus ~ - 0.00002,
  nla ~ 73.046,
  nus ~ 73.054,
  r_pdsus ~ 0.95648,
  r_pdsla ~ 0.95649,
  pmla ~ 0.95628,
  pmus ~ 0.95661,
  r_plas ~ 0.95646,
  r_puss ~ 0.9565,
  pxla ~ 0.95634,
  pxus ~ 0.95656,
  r_play ~ 0.95648,
  r_pusy ~ 0.95649,
  r_sla ~ 109.32,
  r_sus ~ 109.33,
  sla ~ 104.56,
  sus ~ 104.57,
  tla ~ 19.463,
  tus ~ 19.465,
  r_xla ~ 11.926,
  r_xus ~ 11.928,
  xla ~ 11.406,
  xus ~ 11.41,
  xrla ~ 1.0003,
  xrus ~ 0.99971,
  #xrlae ~ 1.0003,
  #xruse ~ 0.99971,
  r_yla ~ 97.392,
  r_yus ~ 97.403,
  yla ~ 93.154,
  yus ~ 93.164,
  ydla ~ 77.851,
  ydus ~ 77.86,
  r_ydla ~ 81.394,
  r_ydus ~ 81.402,
  r_ydlae ~ 81.394,
  r_yduse ~ 81.402,
  rerus~0,
  rerla~0,
  rdus~0.01,
  rdla~0.01,
  pcba~ 6.25,
  pder~ 6.25,
  dxrlae~0,
  dxruse~0,
  tbla~0,
  tbus~0,
  psbrla~0,
  psbrus~0,
  prbla~0,
  prbus~0,
  nwla~0,
  nwus~0,
  rsh_cla~0,
  sh_cabla~0,
  sh_tbla~0,
  sh_govdef~0,
  sh_prbla~0,
  nip~0,
  govdeb~0,
  sh_bbusla_s~0,
  sh_cba~0,
  totla~0
)

cba <- sfcr_baseline(
  equations = eqs, 
  external = external,
  initial=initial,
  periods = 150, 
  
)

#cba

shock1 <- sfcr_shock(
  variables = sfcr_set(
    e_la ~ - 1.08
  ),
  start = 50,
  end = 60
  
)

shock2 <- sfcr_shock(
  variables = sfcr_set(
    e_la ~ - 1.28
  ),
  start = 61,
  end = 150
  
)

cba2 <- sfcr_scenario(
  baseline = cba,
  scenario = list(shock1, shock2),
  periods = 150
)

#cba2

sfcr_dag_blocks(eqs)

sfcr_dag_cycles(eqs)


# Balance-sheet matrix

bs <- sfcr_matrix(
  columns = c("Households CHL", "Firms CHL", "Financial Sector CHL", "Government CHL", "Central bank CHL", "FX", "Households RoW", "Firms RoW", "Financial Sector RoW", "Government RoW", "Central bank RoW", "sum"),
  codes = c("hCHL", "fCHL", "bCHL", "gCHL", "cbCHL", "xr", "hRoW", "fRoW", "bRoW", "gRoW", "cbRoW", "s"),
  r1 = c("Deposits", hCHL = "+D", bCHL = "-D", xr = "xr", hRoW = "+D", bRoW = "-D"),
  r2 = c("Bills CHL", hCHL = "+B_hCHLCHL", bCHL = "+B_bCHLCHL", gCHL = "-B_CHL", cbCHL = "+B_cbCHL", xr = "xr", bRoW = "+B_hRoWCHL"),
  r3 = c("Bills RoW", bCHL = "+B_bCHLRoW", cbCHL = "+B_cbCHLRoW", xr = "xr", hRoW = "+B_hRoWRoW", bRoW = "+B_bRoWRoW", gRoW = "-B_RoW", cbRoW = "+B_cbRoW"),
  r4 = c("Derivatives", hCHL = "+pDER", bCHL = "-pDER", xr = "xr"),
  r5 = c("CLNs", xr = "xr", hRoW = "+pCLN", bRoW = "-pCLN"),
  r6 = c("HPM", hCHL = "+H_h", bCHL = "+H_b", cbCHL = "-H", xr = "xr", hRoW = "+H_h", bRoW = "+H_b", cbRoW = "-H"),
  r7 = c("Balance", hCHL = "-V", bCHL = "0", gCHL = "-NW_g",  cbCHL = "-NW_cb", hRoW = "-V", bRoW = "0", gRoW = "-NW_g",  cbRoW = "-NW_cb")
)


# Transactions-flow matrix (with xr sankey not working, hence I excluded it)
tfm <- sfcr_matrix(
  columns = c("Households CHL", "Firms CHL", "Financial Sector CHL", "Government CHL", "Central bank CHL", "FX", "Households RoW", "Firms RoW", "Financial Sector RoW", "Government RoW", "Central bank RoW"),
  codes = c("hCHL", "fCHL", "bCHL", "gCHL", "cbCHL", "xr", "hRoW", "fRoW", "bRoW", "gRoW", "cbRoW"),
  c("Consumption", hCHL = "-cla", fCHL = "+cla", xr = "xr", hRoW = "-cus", fRoW = "+cus"),
  c("Govt. Expenditures", fCHL = "+gla", gCHL = "-gla", xr = "xr", fRoW = "+gus", gRoW = "-gus"),
  c("Exports", fCHL = "+xla", xr = "xr", fRoW = "-imus"),
  c("Imports", fCHL = "-imla", xr = "xr", fRoW = "+xus"),
  c("Taxes", fCHL = "-tla", gCHL = "+tla", xr = "xr", fRoW = "-tus", gRoW = "+tus"),
  c("Income", hCHL = "+yla", fCHL = "-yla", xr = "xr", hRoW = "+yus", fRoW = "-yus"),
  c("Banks Profits", hCHL = "+fbla", bCHL = "-fbla", xr = "xr", hRoW = "+fbus", bRoW = "-fbus"),
  c("CB Profits", gCHL = "+fcbla", cbCHL = "-fcbla", xr = "xr", gRoW = "+fcbus", cbRoW = "-fcbus"),
  c("Int. Deposits", hCHL = "+rdla[-1]*dla_d[-1]", bCHL = "-rdla[-1]*dla_d[-1]", xr = "xr", hRoW = "+rdus[-1]*dus_d[-1]", bRoW = "-rdus[-1]*dus_d[-1]"),
  c("Int. Bills CHL", hCHL = "+rla[-1]*blala_d[-1]", bCHL = "+rla[-1]*bbla_d[-1]", gCHL = "-rla[-1]*bla_s[-1]", cbCHL = "+rla[-1]*bcbla_d[-1]", xr = "xr", bRoW = "+rla[-1]*bbusla_d[-1]"),
  c("Int. Bills RoW", bCHL = "+rus[-1]*bblaus_d[-1]", cbCHL = "+rus[-1]*bcblaus_d[-1]", xr = "xr", hRoW = "+rus[-1]*busus_d[-1]", bRoW = "+rus[-1]*bbus_d[-1]", gRoW = "-rus[-1]*bus_s[-1]", cbRoW = "+rus[-1]*bcbus_d[-1]"),
  c("Int. Derivatives", hCHL = "+rder[-1]*der_d[-1]", bCHL = "-rder[-1]*der_d[-1]"),
  c("Int. CLNs", hRoW = "+rcba[-1]*cba_d[-1]", bRoW = "-rcba[-1]*cba_d[-1]"),
  c("Ch. Deposits", hCHL = "-(dla_s-dla_s[-1])", bCHL = "+(dla_s-dla_s[-1])", xr = "xr", hRoW = "-(dus_s-dus_s[-1])", bRoW = "+(dus_s-dus_s[-1])"),
  c("Ch. Bills CHL", hCHL = "-(blala_s-blala_s[-1])", bCHL = "-(bbla_s-bbla_s[-1])", gCHL = "+(bla_s-bla_s[-1])", cbCHL = "-(bcbla_s-bcbla_s[-1])", xr = "xr", bRoW = "-(bbusla_s-bbusla_s[-1])"),
  c("Ch. Bills RoW", bCHL = "-(bblaus_s-bblaus_s[-1])", cbCHL = "-(bcblaus_s-bcblaus_s[-1])", xr = "xr", hRoW = "-(busus_s-busus_s[-1])", bRoW = "-(bbus_s-bbus_s[-1])", gRoW = "+(bus_s-bus_s[-1])", cbRoW = "-(bcbus_s-bcbus_s[-1])"),
  c("Ch. Derivatives", hCHL = "-(pder*der_s-pder[-1]*der_s[-1])", bCHL = "+(pder*der_s-pder[-1]*der_s[-1])" ),
  c("Ch. CLNs", hRoW = "-(pcba*cba_s-pcba[-1]*cba_s[-1])", bRoW = "+(pcba*cba_s-pcba[-1]*cba_s[-1])"),
  c("Ch. HPM", hCHL = "-(hla_s-hla_s[-1])", bCHL = "-(hbla_s-hbla_s[-1])", cbCHL = "+(hla_s-hla_s[-1])+(hbla_s-hbla_s[-1])", xr = "xr", hRoW = "-(hus_s-hus_s[-1])", bRoW = "-(hbus_s-hbus_s[-1])", cbRoW = "+(hus_s-hus_s[-1])+(hbus_s-hbus_s[-1])"),
)


p1 <- cba2 %>%
  ggplot( aes(x=period, y=r_yla)) +
  geom_line(aes(colour="CHL real output")) +
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())
p2 <- cba2 %>%
  ggplot( aes(x=period, y=rsh_cla)) +
  geom_line(aes(colour="CHL consumption % GDP")) +
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)
p3 <- cba2 %>%
  ggplot( aes(x=period, y=sh_cabla)) +
  geom_line( aes(color="CA % GDP"))+
  geom_line(aes(x=period, y=sh_tbla, color="TA % GDP"))  +
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))
p4 <- cba2 %>%
  ggplot( aes(x=period, y=sh_cabla)) +
  geom_line(aes(color="CA % GDP")) +
  geom_line(aes(x=period, y=sh_govdef, color="GB % GDP")) +
  geom_line(aes(x=period, y=sh_prbla, color="PB % GDP")) +
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values =c("blue", "red", "green") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))

grid.arrange(p1, p2, p3, p4, nrow = 2)

p5 <- cba2 %>%
  ggplot( aes(x=period, y = rerla)) +
  geom_line(aes(color="RER CHL")) +
  geom_line(aes(y = pxla/pmla, color="TOT CHL"))  +
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))
p6 <- cba2 %>%
  ggplot( aes(x=period, y=cgus)) +
  geom_line(aes(color="Capital gains on CLNs")) +
  geom_line(aes(x=period, y=cgbus, color="Capital gains on CHL bills"))   +
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))
p7 <- cba2 %>%
  ggplot( aes(x=period, y=pcba)) +
  geom_line(aes(color="Price of CLNs"))  +
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) 
p8 <- cba2 %>%
  ggplot( aes(x=period, y=nip)) +
  geom_line(aes(color="Interest payment on CHL debt % GDP"))  +
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) 

grid.arrange(p5, p6, p7, p8, nrow = 2)

p9 <- cba2 %>%
  ggplot( aes(x=period, y = govdeb)) +
  geom_bar(stat = "identity", aes(color="Government Debt")) +
  coord_cartesian(ylim=c(1.53,1.60)) +
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) 
p10 <- cba2 %>%
  ggplot( aes(x=period, y=nwla)) +
  #geom_bar(aes(x=period, y=nwus)) +
  geom_bar(stat = "identity", aes(color="Neth Wealth CHL")) +
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())
p11 <- cba2 %>%
  ggplot( aes(x=period, y=sh_bbusla_s)) +
  geom_bar(stat = "identity", aes(color="External to Total Gov. Debt")) +
  coord_cartesian(ylim=c(0.06,0.092)) +
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) 
p12 <- cba2 %>%
  ggplot( aes(x=period, y=sh_cba)) +
  geom_bar(stat = "identity", aes(color="Share of CLN in RoW wealth")) +
  coord_cartesian(ylim=c(0.150,0.162))+
  theme(axis.title.x=element_blank())+
  theme(axis.title.y=element_blank()) +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) 

grid.arrange(p9, p10, p11, p12, nrow = 2)

# Block and cycles

#sfcr_dag_blocks_plot(eqs, title = NULL, size = 10)


#sfcr_dag_cycles_plot(eqs, title = NULL, size = 10)

#sfcr_sankey(tfm, cba, when = "start")

#Sensitivity analysis: fractional Reserve of RoW intermediaries.
external2 <- sfcr_set(
  alpha_1la ~ 0.9,
  alpha_1us ~ 0.6,
  alpha_2la ~ 0.053333276,
  alpha_2us ~ 0.213326723889398,
  #vertical and horizontal constraints respected
  lambda_10la ~ 0.4,
  lambda_11la ~ 0.5,
  lambda_12la ~ 0.5,
  lambda_13la~ 0.25,
  lambda_14la~ 0.25,
  lambda_20la ~ 0.3,
  lambda_21la ~ 0.5,
  lambda_22la ~ 0.5,
  lambda_23la~ 0.5,
  lambda_40la ~ 0.4,
  lambda_41la ~ 0.25,
  lambda_42la ~ 5,
  lambda_43la~ 0.5,
  lambda_44la~ 0.25,
  lambda_50la ~ 0.1,
  lambda_51la ~ 0.25,
  lambda_52la ~ 5,
  lambda_53la~ 0.25,
  lambda_54la~ 0.5,
  #vertical and horizontal constraints respected
  lambda_24la~ 5,
  lambda_30la ~ 0.25,
  lambda_31la ~ 5,
  lambda_32la ~ 5,
  lambda_33la~ 5,
  lambda_34la~ 5,
  #vertical and horizontal constraints respected
  lambda_10us ~ 0.4,
  lambda_11us ~ 0.5,
  lambda_12us ~ 0.5,
  lambda_13us~ 0.25,
  lambda_14us~ 0.25,
  lambda_20us ~ 0.3,
  lambda_21us ~ 0.5,
  lambda_22us ~ 0.5,
  lambda_23us~ 0.5,
  lambda_40us ~ 0.4,
  lambda_41us ~ 0.25,
  lambda_42us ~ 5,
  lambda_43us~ 0.5,
  lambda_44us~ 0.25,
  lambda_50us ~ 0.1,
  lambda_51us ~ 0.25,
  lambda_52us ~ 5,
  lambda_53us~ 0.25,
  lambda_54us~ 0.5,
  #vertical and horizontal constraints respected
  lambda_24us~ 5,
  lambda_30us ~ 0.25,
  lambda_31us ~ 5,
  lambda_32us ~ 5,
  lambda_33us~ 5,
  lambda_34us~ 5,
  #vertical and horizontal constraints respected
  lambda_60us~0.5,
  lambda_61us~0.5,
  lambda_62us~0.5,
  lambda_70us~0.5,
  lambda_71us~0.5,
  lambda_72us~0.5,
  lambda_60la~0.5,
  lambda_61la~0.5,
  lambda_62la~0.5,
  lambda_70la~0.5,
  lambda_71la~0.5,
  lambda_72la~0.5,
  #vertical and horizontal constraints respected
  e_la ~ - 1.18,
  eta_la ~ 0.7,
  epsilon_la ~ 0.8,
  p_la ~ - 3.015,
  psi_la ~ 0.7,
  pi_la ~  1.2,
  mi_0la ~ - 0.00001,
  chi_0la ~ - 0.00001,
  mi_1la ~ 0.7,
  chi_1la ~ 0.2,
  phila ~ 0.2381,
  phius ~ 0.2381,
  thetala ~ 0.2,
  thetaus ~ 0.2,
  thetabla ~ 0,
  thetabus ~ 0,
  rho_0us~0.05,
  rho_1us~0.9,
  rho_2us~2.7,
  rho_0la~0.03,
  rho_1la~0.9,
  rho_2la~2.7,
  rho_0int~0.9,
  rho_1int~0.9,
  rho_2int~0.9,
  # exogenous variables
  bcblaus_s ~ 0.02031,
  r_gla ~ 16,
  r_gus ~ 16,
  r_prla ~ 1.3333,
  r_prus ~ 1.3333,
  rla ~ 0.03,
  rus ~ 0.03,
  rcba~0.16,
  rder~0.16,
  wla ~ 1,
  wus ~ 1
)

cba9 <- sfcr_baseline(
  equations = eqs, 
  external = external2,
  initial=initial,
  periods = 150, 
  
)

#cba3


cba10 <- sfcr_scenario(
  baseline = cba9,
  scenario = list(shock1, shock2),
  periods = 150
)

#cba4

p1 <- cba10 %>%
  ggplot( aes(x=period, y=r_yla)) +
  geom_line(aes(colour="CHL real output")) +
  theme(axis.title.x=element_blank())+
  labs(y ="ychl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
p2 <- cba10 %>%
  ggplot( aes(x=period, y=rsh_cla)) +
  geom_line(aes(colour="CHL consumption % GDP")) +
  theme(axis.title.x=element_blank())+
  labs(y ="Cchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = "none")
p3 <- cba10 %>%
  ggplot( aes(x=period, y=sh_cabla)) +
  geom_line( aes(color="CA % GDP"))+
  #geom_line(aes(x=period, y=sh_tbla, color="TA % GDP"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="CA") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p4 <- cba10 %>%
  ggplot( aes(x=period, y=sh_govdef)) +
  geom_line(aes(color="CA % GDP")) +
  #geom_line(aes(x=period, y=sh_govdef, color="GB % GDP")) +
  #geom_line(aes(x=period, y=sh_prbla, color="PB % GDP")) +
  theme(axis.title.x=element_blank())+
  labs(y ="GB") +
  scale_color_manual(values =c("blue", "red", "green") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
#grid.arrange(p1, p2, p3, p4, nrow = 2)

p5 <- cba10 %>%
  ggplot( aes(x=period, y = rerla)) +
  geom_line(aes(color="RER CHL")) +
  #geom_line(aes(y = pxla/pmla, color="TOT CHL"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="RERchl") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p6 <- cba10 %>%
  ggplot( aes(x=period, y=cgus)) +
  geom_line(aes(color="Capital gains on CLNs")) +
  #geom_line(aes(x=period, y=cgbus, color="Capital gains on CHL bills"))   +
  theme(axis.title.x=element_blank())+
  labs(y ="GCcln") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p7 <- cba10 %>%
  ggplot( aes(x=period, y=pcba)) +
  geom_line(aes(color="Price of CLNs"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="pCLN") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
p8 <- cba10 %>%
  ggplot( aes(x=period, y=nip)) +
  geom_line(aes(color="Interest payment on CHL debt % GDP"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="nipchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

#grid.arrange(p5, p6, p7, p8, nrow = 2)

p9 <- cba10 %>%
  ggplot( aes(x=period, y = govdeb)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(1.50,1.65)) +
  theme(axis.title.x=element_blank())+
  labs(y ="govdebchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")
p10 <- cba10 %>%
  ggplot( aes(x=period, y=nwla)) +
  #geom_bar(aes(x=period, y=nwus)) +
  geom_bar(stat = "identity", color="red") +
  theme(axis.title.x=element_blank())+
  labs(y ="nwchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  theme(legend.position = "none")
p11 <- cba10 %>%
  ggplot( aes(x=period, y=sh_bbusla_s)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(0.04,0.11)) +
  theme(axis.title.x=element_blank())+
  labs(y ="sh_bbrowchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")
p12 <- cba10 %>%
  ggplot( aes(x=period, y=sh_cba)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(0.147,0.165))+
  theme(axis.title.x=element_blank())+
  labs(y ="sh_cln") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

#grid.arrange(p9, p10, p11, p12, nrow = 2)

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow = 4)

#Sensitivity analysis: fraction of Foreign Bills in RoW Intermediaries’ portfolio
external3 <- sfcr_set(
  alpha_1la ~ 0.9,
  alpha_1us ~ 0.6,
  alpha_2la ~ 0.053333276,
  alpha_2us ~ 0.213326723889398,
  #vertical and horizontal constraints respected
  lambda_10la ~ 0.4,
  lambda_11la ~ 0.5,
  lambda_12la ~ 0.5,
  lambda_13la~ 0.25,
  lambda_14la~ 0.25,
  lambda_20la ~ 0.3,
  lambda_21la ~ 0.5,
  lambda_22la ~ 0.5,
  lambda_23la~ 0.5,
  lambda_40la ~ 0.4,
  lambda_41la ~ 0.25,
  lambda_42la ~ 5,
  lambda_43la~ 0.5,
  lambda_44la~ 0.25,
  lambda_50la ~ 0.1,
  lambda_51la ~ 0.25,
  lambda_52la ~ 5,
  lambda_53la~ 0.25,
  lambda_54la~ 0.5,
  #vertical and horizontal constraints respected
  lambda_24la~ 5,
  lambda_30la ~ 0.25,
  lambda_31la ~ 5,
  lambda_32la ~ 5,
  lambda_33la~ 5,
  lambda_34la~ 5,
  #vertical and horizontal constraints respected
  lambda_10us ~ 0.4,
  lambda_11us ~ 0.5,
  lambda_12us ~ 0.5,
  lambda_13us~ 0.25,
  lambda_14us~ 0.25,
  lambda_20us ~ 0.3,
  lambda_21us ~ 0.5,
  lambda_22us ~ 0.5,
  lambda_23us~ 0.5,
  lambda_40us ~ 0.4,
  lambda_41us ~ 0.25,
  lambda_42us ~ 5,
  lambda_43us~ 0.5,
  lambda_44us~ 0.25,
  lambda_50us ~ 0.1,
  lambda_51us ~ 0.25,
  lambda_52us ~ 5,
  lambda_53us~ 0.25,
  lambda_54us~ 0.5,
  #vertical and horizontal constraints respected
  lambda_24us~ 5,
  lambda_30us ~ 0.25,
  lambda_31us ~ 5,
  lambda_32us ~ 5,
  lambda_33us~ 5,
  lambda_34us~ 5,
  #vertical and horizontal constraints respected
  lambda_60us~0.5,
  lambda_61us~0.5,
  lambda_62us~0.5,
  lambda_70us~0.5,
  lambda_71us~0.5,
  lambda_72us~0.5,
  lambda_60la~0.5,
  lambda_61la~0.5,
  lambda_62la~0.5,
  lambda_70la~0.5,
  lambda_71la~0.5,
  lambda_72la~0.5,
  #vertical and horizontal constraints respected
  e_la ~ - 1.18,
  eta_la ~ 0.7,
  epsilon_la ~ 0.8,
  p_la ~ - 3.015,
  psi_la ~ 0.7,
  pi_la ~  1.2,
  mi_0la ~ - 0.00001,
  chi_0la ~ - 0.00001,
  mi_1la ~ 0.7,
  chi_1la ~ 0.2,
  phila ~ 0.2381,
  phius ~ 0.2381,
  thetala ~ 0.2,
  thetaus ~ 0.2,
  thetabla ~ 0,
  thetabus ~ 0,
  rho_0us~0.03,
  rho_1us~0.9,
  rho_2us~1.35,
  rho_0la~0.03,
  rho_1la~0.9,
  rho_2la~2.7,
  rho_0int~0.9,
  rho_1int~0.9,
  rho_2int~0.9,
  # exogenous variables
  bcblaus_s ~ 0.02031,
  r_gla ~ 16,
  r_gus ~ 16,
  r_prla ~ 1.3333,
  r_prus ~ 1.3333,
  rla ~ 0.03,
  rus ~ 0.03,
  rcba~0.16,
  rder~0.16,
  wla ~ 1,
  wus ~ 1
)

cba11 <- sfcr_baseline(
  equations = eqs, 
  external = external3,
  initial=initial,
  periods = 150, 
  
)

#cba3


cba12 <- sfcr_scenario(
  baseline = cba11,
  scenario = list(shock1, shock2),
  periods = 150
)

#cba4

p1 <- cba12 %>%
  ggplot( aes(x=period, y=r_yla)) +
  geom_line(aes(colour="CHL real output")) +
  theme(axis.title.x=element_blank())+
  labs(y ="ychl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
p2 <- cba12 %>%
  ggplot( aes(x=period, y=rsh_cla)) +
  geom_line(aes(colour="CHL consumption % GDP")) +
  theme(axis.title.x=element_blank())+
  labs(y ="Cchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = "none")
p3 <- cba12 %>%
  ggplot( aes(x=period, y=sh_cabla)) +
  geom_line( aes(color="CA % GDP"))+
  #geom_line(aes(x=period, y=sh_tbla, color="TA % GDP"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="CA") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p4 <- cba12 %>%
  ggplot( aes(x=period, y=sh_govdef)) +
  geom_line(aes(color="CA % GDP")) +
  #geom_line(aes(x=period, y=sh_govdef, color="GB % GDP")) +
  #geom_line(aes(x=period, y=sh_prbla, color="PB % GDP")) +
  theme(axis.title.x=element_blank())+
  labs(y ="GB") +
  scale_color_manual(values =c("blue", "red", "green") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
#grid.arrange(p1, p2, p3, p4, nrow = 2)

p5 <- cba12 %>%
  ggplot( aes(x=period, y = rerla)) +
  geom_line(aes(color="RER CHL")) +
  #geom_line(aes(y = pxla/pmla, color="TOT CHL"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="RERchl") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p6 <- cba12 %>%
  ggplot( aes(x=period, y=cgus)) +
  geom_line(aes(color="Capital gains on CLNs")) +
  #geom_line(aes(x=period, y=cgbus, color="Capital gains on CHL bills"))   +
  theme(axis.title.x=element_blank())+
  labs(y ="GCcln") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p7 <- cba12 %>%
  ggplot( aes(x=period, y=pcba)) +
  geom_line(aes(color="Price of CLNs"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="pCLN") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
p8 <- cba12 %>%
  ggplot( aes(x=period, y=nip)) +
  geom_line(aes(color="Interest payment on CHL debt % GDP"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="nipchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

#grid.arrange(p5, p6, p7, p8, nrow = 2)

p9 <- cba12 %>%
  ggplot( aes(x=period, y = govdeb)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(1.4,1.55)) +
  theme(axis.title.x=element_blank())+
  labs(y ="govdebchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")
p10 <- cba12 %>%
  ggplot( aes(x=period, y=nwla)) +
  #geom_bar(aes(x=period, y=nwus)) +
  geom_bar(stat = "identity", color="red") +
  theme(axis.title.x=element_blank())+
  labs(y ="nwchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  theme(legend.position = "none")
p11 <- cba12 %>%
  ggplot( aes(x=period, y=sh_bbusla_s)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(0.02,0.05)) +
  theme(axis.title.x=element_blank())+
  labs(y ="sh_bbrowchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")
p12 <- cba12 %>%
  ggplot( aes(x=period, y=sh_cba)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(0.147,0.165))+
  theme(axis.title.x=element_blank())+
  labs(y ="sh_cln") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

#grid.arrange(p9, p10, p11, p12, nrow = 2)

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow = 4)

#Sensitivity analysis: spread between CHL deposit rate and Bills.
external4 <- sfcr_set(
  alpha_1la ~ 0.9,
  alpha_1us ~ 0.6,
  alpha_2la ~ 0.053333276,
  alpha_2us ~ 0.213326723889398,
  #vertical and horizontal constraints respected
  lambda_10la ~ 0.4,
  lambda_11la ~ 0.5,
  lambda_12la ~ 0.5,
  lambda_13la~ 0.25,
  lambda_14la~ 0.25,
  lambda_20la ~ 0.3,
  lambda_21la ~ 0.5,
  lambda_22la ~ 0.5,
  lambda_23la~ 0.5,
  lambda_40la ~ 0.4,
  lambda_41la ~ 0.25,
  lambda_42la ~ 5,
  lambda_43la~ 0.5,
  lambda_44la~ 0.25,
  lambda_50la ~ 0.1,
  lambda_51la ~ 0.25,
  lambda_52la ~ 5,
  lambda_53la~ 0.25,
  lambda_54la~ 0.5,
  #vertical and horizontal constraints respected
  lambda_24la~ 5,
  lambda_30la ~ 0.25,
  lambda_31la ~ 5,
  lambda_32la ~ 5,
  lambda_33la~ 5,
  lambda_34la~ 5,
  #vertical and horizontal constraints respected
  lambda_10us ~ 0.4,
  lambda_11us ~ 0.5,
  lambda_12us ~ 0.5,
  lambda_13us~ 0.25,
  lambda_14us~ 0.25,
  lambda_20us ~ 0.3,
  lambda_21us ~ 0.5,
  lambda_22us ~ 0.5,
  lambda_23us~ 0.5,
  lambda_40us ~ 0.4,
  lambda_41us ~ 0.25,
  lambda_42us ~ 5,
  lambda_43us~ 0.5,
  lambda_44us~ 0.25,
  lambda_50us ~ 0.1,
  lambda_51us ~ 0.25,
  lambda_52us ~ 5,
  lambda_53us~ 0.25,
  lambda_54us~ 0.5,
  #vertical and horizontal constraints respected
  lambda_24us~ 5,
  lambda_30us ~ 0.25,
  lambda_31us ~ 5,
  lambda_32us ~ 5,
  lambda_33us~ 5,
  lambda_34us~ 5,
  #vertical and horizontal constraints respected
  lambda_60us~0.5,
  lambda_61us~0.5,
  lambda_62us~0.5,
  lambda_70us~0.5,
  lambda_71us~0.5,
  lambda_72us~0.5,
  lambda_60la~0.5,
  lambda_61la~0.5,
  lambda_62la~0.5,
  lambda_70la~0.5,
  lambda_71la~0.5,
  lambda_72la~0.5,
  #vertical and horizontal constraints respected
  e_la ~ - 1.18,
  eta_la ~ 0.7,
  epsilon_la ~ 0.8,
  p_la ~ - 3.015,
  psi_la ~ 0.7,
  pi_la ~  1.2,
  mi_0la ~ - 0.00001,
  chi_0la ~ - 0.00001,
  mi_1la ~ 0.7,
  chi_1la ~ 0.2,
  phila ~ 0.2381,
  phius ~ 0.2381,
  thetala ~ 0.2,
  thetaus ~ 0.2,
  thetabla ~ 0,
  thetabus ~ 0,
  rho_0us~0.03,
  rho_1us~0.9,
  rho_2us~2.7,
  rho_0la~0.03,
  rho_1la~0.45,
  rho_2la~2.7,
  rho_0int~0.9,
  rho_1int~0.9,
  rho_2int~0.9,
  # exogenous variables
  bcblaus_s ~ 0.02031,
  r_gla ~ 16,
  r_gus ~ 16,
  r_prla ~ 1.3333,
  r_prus ~ 1.3333,
  rla ~ 0.03,
  rus ~ 0.03,
  rcba~0.16,
  rder~0.16,
  wla ~ 1,
  wus ~ 1
)

cba13 <- sfcr_baseline(
  equations = eqs, 
  external = external4,
  initial=initial,
  periods = 150, 
  
)

#cba3


cba14 <- sfcr_scenario(
  baseline = cba13,
  scenario = list(shock1, shock2),
  periods = 150
)

#cba4

p1 <- cba14 %>%
  ggplot( aes(x=period, y=r_yla)) +
  geom_line(aes(colour="CHL real output")) +
  theme(axis.title.x=element_blank())+
  labs(y ="ychl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
p2 <- cba14 %>%
  ggplot( aes(x=period, y=rsh_cla)) +
  geom_line(aes(colour="CHL consumption % GDP")) +
  theme(axis.title.x=element_blank())+
  labs(y ="Cchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = "none")
p3 <- cba14 %>%
  ggplot( aes(x=period, y=sh_cabla)) +
  geom_line( aes(color="CA % GDP"))+
  #geom_line(aes(x=period, y=sh_tbla, color="TA % GDP"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="CA") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p4 <- cba14 %>%
  ggplot( aes(x=period, y=sh_govdef)) +
  geom_line(aes(color="CA % GDP")) +
  #geom_line(aes(x=period, y=sh_govdef, color="GB % GDP")) +
  #geom_line(aes(x=period, y=sh_prbla, color="PB % GDP")) +
  theme(axis.title.x=element_blank())+
  labs(y ="GB") +
  scale_color_manual(values =c("blue", "red", "green") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
#grid.arrange(p1, p2, p3, p4, nrow = 2)

p5 <- cba14 %>%
  ggplot( aes(x=period, y = rerla)) +
  geom_line(aes(color="RER CHL")) +
  #geom_line(aes(y = pxla/pmla, color="TOT CHL"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="RERchl") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p6 <- cba14 %>%
  ggplot( aes(x=period, y=cgus)) +
  geom_line(aes(color="Capital gains on CLNs")) +
  #geom_line(aes(x=period, y=cgbus, color="Capital gains on CHL bills"))   +
  theme(axis.title.x=element_blank())+
  labs(y ="GCcln") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p7 <- cba14 %>%
  ggplot( aes(x=period, y=pcba)) +
  geom_line(aes(color="Price of CLNs"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="pCLN") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
p8 <- cba14 %>%
  ggplot( aes(x=period, y=nip)) +
  geom_line(aes(color="Interest payment on CHL debt % GDP"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="nipchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

#grid.arrange(p5, p6, p7, p8, nrow = 2)

p9 <- cba14 %>%
  ggplot( aes(x=period, y = govdeb)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(1.50,1.65)) +
  theme(axis.title.x=element_blank())+
  labs(y ="govdebchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")
p10 <- cba14 %>%
  ggplot( aes(x=period, y=nwla)) +
  #geom_bar(aes(x=period, y=nwus)) +
  geom_bar(stat = "identity", color="red") +
  theme(axis.title.x=element_blank())+
  labs(y ="nwchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  theme(legend.position = "none")
p11 <- cba14 %>%
  ggplot( aes(x=period, y=sh_bbusla_s)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(0.04,0.11)) +
  theme(axis.title.x=element_blank())+
  labs(y ="sh_bbrowchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")
p12 <- cba14 %>%
  ggplot( aes(x=period, y=sh_cba)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(0.147,0.165))+
  theme(axis.title.x=element_blank())+
  labs(y ="sh_cln") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

#grid.arrange(p9, p10, p11, p12, nrow = 2)

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow = 4)

#Sensitivity analysis:  no ML condition
external1 <- sfcr_set(
  alpha_1la ~ 0.9,
  alpha_1us ~ 0.6,
  alpha_2la ~ 0.053333276,
  alpha_2us ~ 0.213326723889398,
  #vertical and horizontal constraints respected
  lambda_10la ~ 0.4,
  lambda_11la ~ 0.5,
  lambda_12la ~ 0.5,
  lambda_13la~ 0.25,
  lambda_14la~ 0.25,
  lambda_20la ~ 0.3,
  lambda_21la ~ 0.5,
  lambda_22la ~ 0.5,
  lambda_23la~ 0.5,
  lambda_40la ~ 0.4,
  lambda_41la ~ 0.25,
  lambda_42la ~ 5,
  lambda_43la~ 0.5,
  lambda_44la~ 0.25,
  lambda_50la ~ 0.1,
  lambda_51la ~ 0.25,
  lambda_52la ~ 5,
  lambda_53la~ 0.25,
  lambda_54la~ 0.5,
  #vertical and horizontal constraints respected
  lambda_24la~ 5,
  lambda_30la ~ 0.25,
  lambda_31la ~ 5,
  lambda_32la ~ 5,
  lambda_33la~ 5,
  lambda_34la~ 5,
  #vertical and horizontal constraints respected
  lambda_10us ~ 0.4,
  lambda_11us ~ 0.5,
  lambda_12us ~ 0.5,
  lambda_13us~ 0.25,
  lambda_14us~ 0.25,
  lambda_20us ~ 0.3,
  lambda_21us ~ 0.5,
  lambda_22us ~ 0.5,
  lambda_23us~ 0.5,
  lambda_40us ~ 0.4,
  lambda_41us ~ 0.25,
  lambda_42us ~ 5,
  lambda_43us~ 0.5,
  lambda_44us~ 0.25,
  lambda_50us ~ 0.1,
  lambda_51us ~ 0.25,
  lambda_52us ~ 5,
  lambda_53us~ 0.25,
  lambda_54us~ 0.5,
  #vertical and horizontal constraints respected
  lambda_24us~ 5,
  lambda_30us ~ 0.25,
  lambda_31us ~ 5,
  lambda_32us ~ 5,
  lambda_33us~ 5,
  lambda_34us~ 5,
  #vertical and horizontal constraints respected
  lambda_60us~0.5,
  lambda_61us~0.5,
  lambda_62us~0.5,
  lambda_70us~0.5,
  lambda_71us~0.5,
  lambda_72us~0.5,
  lambda_60la~0.5,
  lambda_61la~0.5,
  lambda_62la~0.5,
  lambda_70la~0.5,
  lambda_71la~0.5,
  lambda_72la~0.5,
  #vertical and horizontal constraints respected
  e_la ~ - 1.18,
  eta_la ~ 0.5,
  epsilon_la ~ 0.8,
  p_la ~ - 3.015,
  psi_la ~ 0.5,
  pi_la ~  1.2,
  mi_0la ~ - 0.00001,
  chi_0la ~ - 0.00001,
  mi_1la ~ 0.7,
  chi_1la ~ 0.2,
  phila ~ 0.2381,
  phius ~ 0.2381,
  thetala ~ 0.2,
  thetaus ~ 0.2,
  thetabla ~ 0,
  thetabus ~ 0,
  rho_0us~0.03,
  rho_1us~0.9,
  rho_2us~2.7,
  rho_0la~0.03,
  rho_1la~0.9,
  rho_2la~2.7,
  rho_0int~0.9,
  rho_1int~0.9,
  rho_2int~0.9,
  # exogenous variables
  bcblaus_s ~ 0.02031,
  r_gla ~ 16,
  r_gus ~ 16,
  r_prla ~ 1.3333,
  r_prus ~ 1.3333,
  rla ~ 0.03,
  rus ~ 0.03,
  rcba~0.16,
  rder~0.16,
  wla ~ 1,
  wus ~ 1
)

cba3 <- sfcr_baseline(
  equations = eqs, 
  external = external1,
  initial=initial,
  periods = 150, 
  
)

#cba3


cba4 <- sfcr_scenario(
  baseline = cba3,
  scenario = list(shock1, shock2),
  periods = 150
)

#cba4

p1 <- cba4 %>%
  ggplot( aes(x=period, y=r_yla)) +
  geom_line(aes(colour="CHL real output")) +
  theme(axis.title.x=element_blank())+
  labs(y ="ychl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
p2 <- cba4 %>%
  ggplot( aes(x=period, y=rsh_cla)) +
  geom_line(aes(colour="CHL consumption % GDP")) +
  theme(axis.title.x=element_blank())+
  labs(y ="Cchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = "none")
p3 <- cba4 %>%
  ggplot( aes(x=period, y=sh_cabla)) +
  geom_line( aes(color="CA % GDP"))+
  #geom_line(aes(x=period, y=sh_tbla, color="TA % GDP"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="CA") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p4 <- cba4 %>%
  ggplot( aes(x=period, y=sh_govdef)) +
  geom_line(aes(color="CA % GDP")) +
  #geom_line(aes(x=period, y=sh_govdef, color="GB % GDP")) +
  #geom_line(aes(x=period, y=sh_prbla, color="PB % GDP")) +
  theme(axis.title.x=element_blank())+
  labs(y ="GB") +
  scale_color_manual(values =c("blue", "red", "green") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
#grid.arrange(p1, p2, p3, p4, nrow = 2)

p5 <- cba4 %>%
  ggplot( aes(x=period, y = rerla)) +
  geom_line(aes(color="RER CHL")) +
  #geom_line(aes(y = pxla/pmla, color="TOT CHL"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="RERchl") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p6 <- cba4 %>%
  ggplot( aes(x=period, y=cgus)) +
  geom_line(aes(color="Capital gains on CLNs")) +
  #geom_line(aes(x=period, y=cgbus, color="Capital gains on CHL bills"))   +
  theme(axis.title.x=element_blank())+
  labs(y ="GCcln") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p7 <- cba4 %>%
  ggplot( aes(x=period, y=pcba)) +
  geom_line(aes(color="Price of CLNs"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="pCLN") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
p8 <- cba4 %>%
  ggplot( aes(x=period, y=nip)) +
  geom_line(aes(color="Interest payment on CHL debt % GDP"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="nipchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

#grid.arrange(p5, p6, p7, p8, nrow = 2)

p9 <- cba4 %>%
  ggplot( aes(x=period, y = govdeb)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(1.50,1.65)) +
  theme(axis.title.x=element_blank())+
  labs(y ="govdebchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")
p10 <- cba4 %>%
  ggplot( aes(x=period, y=nwla)) +
  #geom_bar(aes(x=period, y=nwus)) +
  geom_bar(stat = "identity", color="red") +
  theme(axis.title.x=element_blank())+
  labs(y ="nwchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  theme(legend.position = "none")
p11 <- cba4 %>%
  ggplot( aes(x=period, y=sh_bbusla_s)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(0.04,0.11)) +
  theme(axis.title.x=element_blank())+
  labs(y ="sh_bbrowchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")
p12 <- cba4 %>%
  ggplot( aes(x=period, y=sh_cba)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(0.147,0.165))+
  theme(axis.title.x=element_blank())+
  labs(y ="sh_cln") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

#grid.arrange(p9, p10, p11, p12, nrow = 2)

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow = 4)

#Sensitivity analysis: different demand for CHL assets for intermediaries
eqs1 <- sfcr_set(
  #model creation,
  #(1) arterial flows,
  #income,
  #yus~ cus+ gus +xus- imus,
  #yla ~ cla+ gla+ xla -imla,
  #disposable income,
  ydus~(yus + fbus +rus[-1]*busus_d[-1]+ rdus[-1]*dus_d[-1] +cba_d[-1]+cgus)*(1-thetaus),
  ydla~(yla + fbla + rla[-1]*blala_d[-1] + rdla[-1]*dla_d[-1]+der_d[-1]+cgla)*(1-thetala), 
  #disposable income (hs),
  #ydus_hs~ydus+ d(xrla)*busla_s[-1],
  #ydla_hs~ydla+ d(xrus)*blaus_s[-1],
  #wealth,
  vus ~ vus[-1]+  ydus-cus,
  vla ~ vla[-1]+ ydla-cla,
  #capital gains,
  cgus ~ (pcba-pcba[-1])*cba_s[-1],
  cgla ~ (pder-pder[-1])*der_s[-1],
  #tax,
  tla ~ thetala*(yla+ fbla+rla[-1]*blala_d[-1]+ rdla[-1]*dla_d[-1]+der_d[-1]),
  tus ~ thetaus*(yus+fbus +rus[-1]*busus_d[-1]+rdus[-1]*dus_d[-1] + cba_d[-1]) ,
  #cb profits,
  fcbla ~ rla[-1]*bcbla_d[-1] + rus[-1]*bcblaus_s[-1]*xrus,
  fcbus ~ rus[-1]*bcbus_d[-1], 
  #government budget constraint,
  bla_s ~ bla_s[-1]+gla-tla+ rla[-1]*bla_s[-1] - fcbla, 
  bus_s ~ bus_s[-1]+gus-tus+ rus[-1]*bus_s[-1] - fcbus,
  #current and capital account,
  cabla~xla - imla - rla[-1]*bbusla_s[-1] + rus[-1]*bblaus_s[-1]*xrus + rus[-1]*bcblaus_s[-1]*xrus ,
  kabla ~ (bbusla_s-bbusla_s[-1]) - (bcblaus_s-bcblaus_s[-1])*xrus- (bblaus_s-bblaus_s[-1])*xrus,
  cabus~xus - imus + rla[-1]*bbusla_s[-1]*xrla - rus[-1]*bblaus_s[-1] - rus[-1]*bcblaus_s[-1],
  #kabus ~ -d(bbla_s)*xrla + d(blaus_s) + d(bcblaus_s),
  #(2) trade,
  #export and import,
  pmla~ exp(mi_0la+  mi_1la *log(r_pusy)  + (1- mi_1la)*log(r_play) - mi_1la *log(xrla)),
  pxla ~exp(chi_0la+  chi_1la *log(r_pusy)  + (1- chi_1la)*log(r_play)- chi_1la *log(xrla)),
  pxus ~ pmla*xrla,  
  pmus ~ pxla*xrla,
  r_xla~ exp(e_la - eta_la *log(pmus/r_pusy) + epsilon_la*log(r_yus)),
  r_imla ~ exp(p_la - psi_la *log(pmla[-1]/r_play[-1]) + pi_la*log(r_yla)),
  r_xus  ~r_imla,
  r_imus~r_xla,
  xla ~ r_xla*pxla,
  xus ~ r_xus*pxus,
  imla ~ r_imla*pmla,  
  imus ~ r_imus*pmus,
  #(3) income and expenditure,
  #real disposable income is of the haig-simon type,
  r_vla ~ vla/r_pdsla, 
  r_vus~ vus/r_pdsus,
  r_ydla ~ ydla/r_pdsla - r_vla[-1]* (r_pdsla-r_pdsla[-1])/r_pdsla ,
  r_ydus ~ ydus/r_pdsus - r_vus[-1]* (r_pdsus-r_pdsus[-1])/r_pdsus ,
  r_cla ~ alpha_1la*r_ydlae + alpha_2la*r_vla[-1],
  r_cus ~alpha_1us*r_yduse + alpha_2us*r_vus[-1],
  r_ydlae~ (r_ydla + r_ydla[-1])/2,
  r_yduse~ (r_ydus + r_ydus[-1])/2,
  r_sla ~r_cla+r_gla+r_xla ,
  r_sus ~r_cus+r_gus+r_xus ,
  sla ~r_sla*r_plas,
  sus ~r_sus*r_puss,
  r_plas ~((1+phila)*(wla*nla+imla)) /r_sla,
  r_puss ~((1+phius)*(wus*nus+imus)) / r_sus,
  r_pdsla ~(sla-xla) / (r_sla-r_xla),
  r_pdsus ~(sus-xus) / (r_sus-r_xus),
  dsla ~sla-xla ,
  dsus ~sus-xus ,
  r_dsla ~r_cla+r_gla ,
  r_dsus ~r_cus+r_gus  ,
  yla~sla-imla,
  yus ~sus-imus,
  r_yla ~r_sla-r_imla ,
  r_yus ~r_sus-r_imus ,
  r_play ~ yla /r_yla ,
  r_pusy ~yus/ r_yus,
  cla ~r_cla*r_pdsla ,
  cus ~r_cus*r_pdsus  ,
  gla ~r_gla*r_pdsla ,
  gus ~r_gus*r_pdsus ,
  nla ~ r_yla/ r_prla , 
  nus ~ r_yus/ r_prus ,
  #(4) financial intermediaries,
  #us financial sector,
  bbus_d~dus_s-bbusla_d-hbus_d+pcba*cba_s,
  #d(vbus)~fbus,
  bbusla_d~rho_2us*cba_s+(1/5)*cgbus,
  hbus_d~rho_0us*dus_s,
  rdus~ rdus[-1] + rho_1us*(rus-rus[-1]),
  cgbus ~ (xrla-xrla[-1])*bbusla_s[-1],
  pcba~(1/rcba)+log(pmus),
  #d(rcba)~rho_0int*,
  #d(rlus)~rho_2us*d(rus),
  fbus~rus[-1]*bbus_d[-1]-rdus*dus_s[-1]+rla[-1]*bbusla_d[-1]*xrla +cgbus-cba_d[-1]-cgus,
  #la financial sector,
  bbla_d~dla_s-bblaus_d-hbla_d+pder*der_s,
  #d(vbla)~fbla,
  bblaus_d~rho_2la*der_s+(1/5)*cgbla,
  hbla_d~rho_0la*dla_s,
  rdla ~ rdla[-1]+ rho_1la*(rla-rla[-1]),
  cgbla ~ (xrus-xrus[-1])*bblaus_s[-1],
  pder~1/rder,
  #d(rder)~-rho_2int*(rla-rus),
  #d(rlla)~rho_2la*d(rla),
  fbla~rla[-1]*bbla_d[-1]-rdla*dla_s[-1]+rus[-1]*bblaus_d[-1]*xrus +cgbla-der_d[-1]-cgla,
  #intermediaries,
  #d(vint)~rus[-1]*bintus_d[-1]+rla[-1]*bintla_d[-1]*xrla +cgint-rcba[-1]*cba_d[-1],
  #cba_s~cba_d,
  #cgint~d(xrla)*bintla_s[-1],
  #pcba~1/rcba,
  #d(rcba)~rho_0int*d(pxla)-rho_1int*(rla-rus),
  #fint~rus[-1]*bintus_d[-1]+rla[-1]*bintla_d[-1]*xrla +cgint-rcba[-1]*cba_d[-1],
  #bintla_d~cba_s-bintus_d,
  #(5) assets demand,
  #asset demand for la resident,
  blala_d~vla*(lambda_10la+lambda_11la*rla-lambda_13la*rdla-lambda_14la*(rder+dxrlae)),
  dla_d~vla*(lambda_40la-lambda_41la*rla+lambda_43la*rdla-lambda_44la*(rder+dxrlae)),
  der_d~(vla/pder)*(lambda_50la-lambda_51la*rla-lambda_53la*rdla+lambda_54la*(rder+dxrlae)),
  #blaus_d~vla*(lambda_20la-lambda_21la*rla-lambda_23la*rdla+lambda_22la*(rus+dxruse)),
  #hla_d/ vla~lambda_30la - lambda_31la*(rus+dxruse) -lambda_32la*rla,
  #asset demand for us resident,
  busus_d~vus*(lambda_10us +lambda_11us*rus-lambda_13us*rdus-lambda_14us*(rcba+dxruse)),
  dus_d~vus*(lambda_40us-lambda_41us*rus+lambda_43us*rdus-lambda_44us*(rcba+dxruse)),
  cba_d~(vus/pcba)*(lambda_50us-lambda_51us*rus-lambda_53us*rdus+lambda_54us*(rcba+dxruse)),
  #busla_d~vus*(lambda_20us -lambda_21us*rus+lambda_22us*(rla+dxrlae)),
  #hus_d/ vus~lambda_30us -lambda_31us*rus -lambda_32us* (rla+dxrlae),
  #asset demand for intermediaries,
  #bintus_d~vint*(lambda_10int + lambda_11int*rus - lambda_12int*(rla+dxrlae)),
  #bintla_d~vint*(lambda_20int - lambda_21int*rus + lambda_22int*(rla+dxrlae)),
  #expected change in the exchange rate (expectations to update),
  dxrlae~d(pder)/ pder,
  dxruse~d(pcba)/ pcba,
  #(6) assets supply,
  #demand for cash,
  hus_d~ vus - busus_d- dus_d- pcba*cba_d ,
  hla_d~ vla - blala_d-dla_d -pder*der_d,
  #cb demand for b, h; d and cba,
  hus_s~ hus_d,
  hla_s~ hla_d,
  dus_s~dus_d,
  dla_s~dla_d,
  busus_s~ busus_d,
  blala_s ~blala_d,
  cba_s~cba_d,
  der_s~der_d,
  bbus_s~ bbus_d,
  bbla_s ~bbla_d,
  #aus_s ~aus_d,
  #ala_s ~ala_d,
  #bintus_s~bintus_d,
  #lus_s~lus_d,
  #lla_s~lla_d,
  hbus_s~hbus_d,
  hbla_s~hbla_d,
  #supply of domestic t bills to cb,
  bcbus_s~ bcbus_d,
  bcbla_s~ bcbla_d ,
  bcbus_d ~ bcbus_d[-1]+ (hus_s-hus_s[-1])+(hbus_s-hbus_s[-1]),
  bcbla_d ~ bcbla_d[-1]+ (hla_s-hla_s[-1])+(hbla_s-hbla_s[-1])-(bcblaus_s-bcblaus_s[-1])*xrus,
  #supply of assets abroad ,
  #busla_s~ bla_s- blala_s- bcbla_s,
  #exchange rate,
  #xrus~ blaus_d /blaus_s,
  xrla~ 1/xrus,
  bblaus_s~bblaus_d*xrla,
  bcblaus_d~bcblaus_s*xrus,
  xrus~ (bbusla_s)/(bbusla_d),
  bbusla_s~ bla_s - blala_s- bcbla_s -  bbla_s,
  #final equation (not possible to appear twice),
  #blaus_s~ blaus_d /xrus,
  rerla~(r_play/r_pusy)*(xrla),
  rerus~1/rerla,
  tbla~xla - imla,
  tbus~xus - imus,
  psbrla~d(bla_s),
  psbrus~d(bus_s),
  prbla~cabla+d(bla_s),
  prbus~cabus+d(bus_s),
  nwla~((vla)-(bla_s) +(bbla_d-dla_s+bblaus_d+hbla_d-pder*der_s)+(bcblaus_d+bcbla_d-hla_s)),
  nwus~((vus)-(bus_s) +(bbus_d-dus_s+bbusla_d+hbus_d-pcba*cba_s)+(bcbus_d-hus_s)),
  rsh_cla~r_cla/yla,
  sh_cabla~cabla/yla,
  sh_tbla~tbla/yla,
  sh_govdef~-psbrla/yla,
  sh_prbla~prbla/yla,
  nip~rla[-1]*bla_s[-1]/yla,
  govdeb~bla_s/yla,
  sh_bbusla_s~bbusla_s/bla_s,
  sh_cba~pcba*cba_d/r_vus,
  totla~pxla/pmla
)

cba5 <- sfcr_baseline(
  equations = eqs1, 
  external = external,
  initial=initial,
  periods = 150, 
  
)

#cba5


cba6 <- sfcr_scenario(
  baseline = cba5,
  scenario = list(shock1, shock2),
  periods = 150
)

#cba6

p1 <- cba6 %>%
  ggplot( aes(x=period, y=r_yla)) +
  geom_line(aes(colour="CHL real output")) +
  theme(axis.title.x=element_blank())+
  labs(y ="ychl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
p2 <- cba6 %>%
  ggplot( aes(x=period, y=rsh_cla)) +
  geom_line(aes(colour="CHL consumption % GDP")) +
  theme(axis.title.x=element_blank())+
  labs(y ="Cchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = "none")
p3 <- cba6 %>%
  ggplot( aes(x=period, y=sh_cabla)) +
  geom_line( aes(color="CA % GDP"))+
  #geom_line(aes(x=period, y=sh_tbla, color="TA % GDP"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="CA") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p4 <- cba6 %>%
  ggplot( aes(x=period, y=sh_govdef)) +
  geom_line(aes(color="CA % GDP")) +
  #geom_line(aes(x=period, y=sh_govdef, color="GB % GDP")) +
  #geom_line(aes(x=period, y=sh_prbla, color="PB % GDP")) +
  theme(axis.title.x=element_blank())+
  labs(y ="GB") +
  scale_color_manual(values =c("blue", "red", "green") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
#grid.arrange(p1, p2, p3, p4, nrow = 2)

p5 <- cba6 %>%
  ggplot( aes(x=period, y = rerla)) +
  geom_line(aes(color="RER CHL")) +
  #geom_line(aes(y = pxla/pmla, color="TOT CHL"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="RERchl") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p6 <- cba6 %>%
  ggplot( aes(x=period, y=cgus)) +
  geom_line(aes(color="Capital gains on CLNs")) +
  #geom_line(aes(x=period, y=cgbus, color="Capital gains on CHL bills"))   +
  theme(axis.title.x=element_blank())+
  labs(y ="GCcln") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p7 <- cba6 %>%
  ggplot( aes(x=period, y=pcba)) +
  geom_line(aes(color="Price of CLNs"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="pCLN") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
p8 <- cba6 %>%
  ggplot( aes(x=period, y=nip)) +
  geom_line(aes(color="Interest payment on CHL debt % GDP"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="nipchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

#grid.arrange(p5, p6, p7, p8, nrow = 2)

p9 <- cba6 %>%
  ggplot( aes(x=period, y = govdeb)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(1.53,1.60)) +
  theme(axis.title.x=element_blank())+
  labs(y ="govdebchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")
p10 <- cba6 %>%
  ggplot( aes(x=period, y=nwla)) +
  #geom_bar(aes(x=period, y=nwus)) +
  geom_bar(stat = "identity", color="red") +
  theme(axis.title.x=element_blank())+
  labs(y ="nwchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  theme(legend.position = "none")
p11 <- cba6 %>%
  ggplot( aes(x=period, y=sh_bbusla_s)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(0.06,0.092)) +
  theme(axis.title.x=element_blank())+
  labs(y ="sh_bbrowchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")
p12 <- cba6 %>%
  ggplot( aes(x=period, y=sh_cba)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(0.150,0.162))+
  theme(axis.title.x=element_blank())+
  labs(y ="sh_cln") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

#grid.arrange(p9, p10, p11, p12, nrow = 2)

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow = 4)

#Sensitivity analysis: simmetric shock


cba7 <- sfcr_baseline(
  equations = eqs, 
  external = external,
  initial=initial,
  periods = 150, 
  
)

#cba5


cba8 <- sfcr_scenario(
  baseline = cba7,
  scenario = shock1,
  periods = 150
)

#cba6

p1 <- cba8 %>%
  ggplot( aes(x=period, y=r_yla)) +
  geom_line(aes(colour="CHL real output")) +
  theme(axis.title.x=element_blank())+
  labs(y ="ychl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
p2 <- cba8 %>%
  ggplot( aes(x=period, y=rsh_cla)) +
  geom_line(aes(colour="CHL consumption % GDP")) +
  theme(axis.title.x=element_blank())+
  labs(y ="Cchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = "none")
p3 <- cba8 %>%
  ggplot( aes(x=period, y=sh_cabla)) +
  geom_line( aes(color="CA % GDP"))+
  #geom_line(aes(x=period, y=sh_tbla, color="TA % GDP"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="CA") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p4 <- cba8 %>%
  ggplot( aes(x=period, y=sh_govdef)) +
  geom_line(aes(color="CA % GDP")) +
  #geom_line(aes(x=period, y=sh_govdef, color="GB % GDP")) +
  #geom_line(aes(x=period, y=sh_prbla, color="PB % GDP")) +
  theme(axis.title.x=element_blank())+
  labs(y ="GB") +
  scale_color_manual(values =c("blue", "red", "green") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent)+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
#grid.arrange(p1, p2, p3, p4, nrow = 2)

p5 <- cba8 %>%
  ggplot( aes(x=period, y = rerla)) +
  geom_line(aes(color="RER CHL")) +
  #geom_line(aes(y = pxla/pmla, color="TOT CHL"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="RERchl") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p6 <- cba8 %>%
  ggplot( aes(x=period, y=cgus)) +
  geom_line(aes(color="Capital gains on CLNs")) +
  #geom_line(aes(x=period, y=cgbus, color="Capital gains on CHL bills"))   +
  theme(axis.title.x=element_blank())+
  labs(y ="GCcln") +
  scale_color_manual(values =c("blue", "red") )+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  guides(colour = guide_legend(nrow = 2))+
  theme(legend.position = "none")
p7 <- cba8 %>%
  ggplot( aes(x=period, y=pcba)) +
  geom_line(aes(color="Price of CLNs"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="pCLN") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank()) +
  theme(legend.position = "none")
p8 <- cba8 %>%
  ggplot( aes(x=period, y=nip)) +
  geom_line(aes(color="Interest payment on CHL debt % GDP"))  +
  theme(axis.title.x=element_blank())+
  labs(y ="nipchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

#grid.arrange(p5, p6, p7, p8, nrow = 2)

p9 <- cba8 %>%
  ggplot( aes(x=period, y = govdeb)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(1.53,1.60)) +
  theme(axis.title.x=element_blank())+
  labs(y ="govdebchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")
p10 <- cba8 %>%
  ggplot( aes(x=period, y=nwla)) +
  #geom_bar(aes(x=period, y=nwus)) +
  geom_bar(stat = "identity", color="red") +
  theme(axis.title.x=element_blank())+
  labs(y ="nwchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  theme(legend.position = "none")
p11 <- cba8 %>%
  ggplot( aes(x=period, y=sh_bbusla_s)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(0.06,0.092)) +
  theme(axis.title.x=element_blank())+
  labs(y ="sh_bbrowchl") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")
p12 <- cba8 %>%
  ggplot( aes(x=period, y=sh_cba)) +
  geom_bar(stat = "identity", color="blue") +
  coord_cartesian(ylim=c(0.150,0.162))+
  theme(axis.title.x=element_blank())+
  labs(y ="sh_cln") +
  scale_color_manual(values = "blue")+ 
  theme(legend.position="bottom") + 
  theme(legend.title = element_blank())+
  scale_y_continuous(labels = scales::percent) +
  theme(legend.position = "none")

#grid.arrange(p9, p10, p11, p12, nrow = 2)

grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow = 4)

