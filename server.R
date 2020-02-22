#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(mrgsolve)
library(deSolve)
library(PKNCA)

code_pk <- '
        $PARAM KA = 0.5, CL = 5, V = 100, Fr = 0.5

        $MAIN
        double CLi = CL*exp(ETA_CL);
        double Vi = V*exp(ETA_V);
        double ke = CLi/Vi;

        $INIT GUT = 0, CENT = 0
        
        $ODE 
        double PKCONC = CENT/Vi*1000;
        dxdt_GUT = -KA*GUT;
        dxdt_CENT = Fr*KA*GUT - ke*CENT;

        $OMEGA @name IIV @labels ETA_CL ETA_V
        0 0
        
        $SIGMA 0
        
        $CAPTURE @annotated
        PKCONC: Concentration (ng/mL)
    '
mod_pk <- mcode("pk_model", code_pk)


code_tumor <- '
        $PARAM KA = 0.5, CL = 5, V = 100, Fr = 0.5, LAMDA0 = 0.005, EMAX = 0.005, EC50 = 50, W0 = 100, FRAC = 0

        $MAIN
        double CLi = CL*exp(ETA_CL);
        double Vi = V*exp(ETA_V);
        double ke = CLi/Vi;
        double EMAXi = EMAX*exp(ETA_EMAX);
        double EC50i = EC50*exp(ETA_EC50);
        
        
        $INIT GUT = 0, CENT = 0, TUMOR1 = 100, TUMOR2 = 0
        
        $ODE 
        dxdt_GUT = -KA*GUT;
        dxdt_CENT = Fr*KA*GUT - ke*CENT;
        
        double PKCONC = CENT/Vi*1000;
        double EFFECT = EMAXi*PKCONC/(EC50i + PKCONC);
        
        dxdt_TUMOR1 = LAMDA0*TUMOR1 - EFFECT*TUMOR1;
        dxdt_TUMOR2 = LAMDA0*TUMOR2;
        
        
        $OMEGA @name IIV @labels ETA_CL ETA_V ETA_EMAX ETA_EC50
        0 0 0 0
        
        $SIGMA 0
        
        $TABLE
        capture TUMOR = TUMOR1 + TUMOR2;

        $CAPTURE PKCONC TUMOR
    '
mod_tumor <- mcode("tumor_growth_inhibition", code_tumor)



code_biomarker <- '
        $PARAM KA = 0.5, CL = 5, V = 100, Fr = 0.5, EMAX = 1, EC50 = 5, GAMMA = 1.5, RECL = 0, REC50 = 0
        
        $MAIN
        double CLi = CL*exp(ETA_CL);
        double Vi = V*exp(ETA_V);
        double ke = CLi/Vi;
        double EMAXi = EMAX*exp(ETA_EMAX);
        double EC50i = EC50*exp(ETA_EC50);
        
        
        $INIT GUT = 0, CENT = 0
        
        $ODE 
        dxdt_GUT = -KA*GUT;
        dxdt_CENT = Fr*KA*GUT - ke*CENT;
        double PKCONC = CENT/Vi*1000;

        $OMEGA @name IIV @labels ETA_CL ETA_V ETA_EMAX ETA_EC50
        0 0 0 0
        
        $SIGMA 0
        
        $TABLE
        capture EFFECT = -EMAXi*pow(PKCONC, GAMMA) / (EC50i + pow(PKCONC, GAMMA));

        $CAPTURE PKCONC EFFECT
    '
mod_biomarker <- mcode("efficacy_biomarker", code_biomarker)



code_anc <- '
        $PARAM KA = 0.5, CL = 5, V = 100, Fr = 0.5, SLOPE = 1.5, CIRC0 = 5, GAMMA = 0.25, MTT = 125

        $MAIN
        double CLi = CL*exp(ETA_CL);
        double Vi = V*exp(ETA_V);
        double ke = CLi/Vi;
        double KTR = 4/MTT;
        double KPROL = KTR;
        double KCIRC = KTR;


        $INIT GUT = 0, CENT = 0, PROL = 5, TRANSIT1 = 5, TRANSIT2 = 5, TRANSIT3 = 5, CIRC = 5

        $ODE
        dxdt_GUT = -KA*GUT;
        dxdt_CENT = Fr*KA*GUT - ke*CENT;

        double PKCONC = CENT/Vi*1000;
        double EFFECT = SLOPE*PKCONC;

        dxdt_PROL = KPROL*PROL*(1-EFFECT)*pow(CIRC0/CIRC, GAMMA) - KTR*PROL;
        
        dxdt_TRANSIT1 = KTR*PROL - KTR*TRANSIT1;
        dxdt_TRANSIT2 = KTR*TRANSIT1 - KTR*TRANSIT2;
        dxdt_TRANSIT3 = KTR*TRANSIT2 - KTR*TRANSIT3;
        
        dxdt_CIRC = KTR*TRANSIT3 - KCIRC*CIRC;


        $OMEGA @name IIV @labels ETA_CL ETA_V
        0 0

        $SIGMA 0

        $TABLE
        capture ANC = CIRC;

        $CAPTURE PKCONC ANC
    '
mod_anc <- mcode("safety_biomarker_anc", code_anc)


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
    
    
    dat_pk <- reactive({
        
        cyc <- input$cycle_pk
        
        amt1 <- input$amt_pk1
        tau1 <- input$interval_pk1
        n1 <- input$n_pk1
        
        amt2 <- input$amt_pk2
        tau2 <- input$interval_pk2
        n2 <- input$n_pk2
        
        cl <-input$cl_pk
        v <- input$v_pk
        ka <- input$ka_pk
        fr <- input$f_pk
        
        dos1 <- expand.ev(cmt = 1,
                          time = c(0,
                                   c(1:cyc) * (tau1*n1 + tau2*n2)),
                          amt = amt1,
                          ii = tau1,
                          addl = n1 - 1,
                          CL = cl,
                          V = v, 
                          KA = ka,
                          Fr = fr,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos2 <- expand.ev(cmt = 1,
                          time = c(c(tau1*n1 + tau2*n2, 
                                     c(1:(cyc)) * (tau1*n1 + tau2*n2)) - tau2*n2)[-1],
                          amt = amt2,
                          ii = tau2,
                          addl = n2 - 1,
                          CL = cl,
                          V = v, 
                          KA = ka,
                          Fr = fr,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos <- expand.ev(cmt = 1,
                         time = 0,
                         amt = 0,
                         ii = 24,
                         addl = 100,
                         CL = cl,
                         V = v, 
                         KA = ka,
                         Fr = fr,
                         ETA_CL = rnorm(1, 0, 0), 
                         ETA_V = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 2)
        
        # union_all(dos1, dos2) %>% arrange(ID, time, addl, ii)
        union_all(dos1, dos2) %>%
            # union_all(dos) %>% 
            arrange(ID, time, addl, ii)
    })
    
    output$distPlotPK <- renderPlot({
        cyc <- input$cycle_pk
        
        amt1 <- input$amt_pk1
        tau1 <- input$interval_pk1
        n1 <- input$n_pk1
        
        amt2 <- input$amt_pk2
        tau2 <- input$interval_pk2
        n2 <- input$n_pk2
        
        mod_pk %>%
            data_set(dat_pk()) %>% 
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>% 
            plot(PKCONC~time)
    })
    
    output$halflife <- renderText({
        # browser()
        round(0.693/(input$cl_pk/input$v_pk), digits = 1)
    })
    
    
    reaccmax <- reactive({
        # PKNCA packages code as followed
        # (since ncappc packages can't be used in linux sys)
        
        cyc <- input$cycle_pk
        
        amt1 <- input$amt_pk1
        tau1 <- input$interval_pk1
        n1 <- input$n_pk1
        
        amt2 <- input$amt_pk2
        tau2 <- input$interval_pk2
        n2 <- input$n_pk2
        
        out <- mod_pk %>%
            data_set(dat_pk()) %>% 
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05)
        
        pk <- out@data %>%
            select(time, PKCONC) %>% 
            mutate(subject = 1, time = time, conc = PKCONC) %>%
            distinct() %>%
            as.data.frame()
        
        # browser()
        
        
        nca <- PKNCAconc(pk, conc~time|subject) %>%
            PKNCAdata(.,
                      intervals=data.frame(
                          start= cyc*(tau1*(n1) + tau2*(n2)) - tau1 - tau2*n2,
                          end= cyc*(tau1*(n1) + tau2*(n2)) - tau2*n2,
                          cmax=TRUE,
                          tmax=TRUE,
                          # aucinf.obs=TRUE,
                          auclast=TRUE
                      )
            ) %>%
            pk.nca() %>%
            as.data.frame(.$result)
        
        nca %>% filter(PPTESTCD == "cmax") %>% .$PPORRES %>% 
                round(digits = 1)
    })
    
    output$pkcmax1 <- renderText({
        # browser()
        reaccmax()
    })
    
    reacauc <- reactive({
        # PKNCA packages code as followed
        # (since ncappc packages can't be used in linux sys)
        
        cyc <- input$cycle_pk
        
        amt1 <- input$amt_pk1
        tau1 <- input$interval_pk1
        n1 <- input$n_pk1
        
        amt2 <- input$amt_pk2
        tau2 <- input$interval_pk2
        n2 <- input$n_pk2
        
        out <- mod_pk %>%
            data_set(dat_pk()) %>% 
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05)
        
        pk <- out@data %>%
            select(time, PKCONC) %>% 
            mutate(subject = 1, time = time, conc = PKCONC) %>%
            distinct() %>%
            as.data.frame()
        
        # browser()
        
        
        nca <- PKNCAconc(pk, conc~time|subject) %>%
            PKNCAdata(.,
                      intervals=data.frame(
                          start= cyc*(tau1*(n1) + tau2*(n2)) - tau1 - tau2*n2,
                          end= cyc*(tau1*(n1) + tau2*(n2)) - tau2*n2,
                          cmax=TRUE,
                          tmax=TRUE,
                          # aucinf.obs=TRUE,
                          auclast=TRUE
                      )
            ) %>%
            pk.nca() %>%
            as.data.frame(.$result)
        
        nca %>% filter(PPTESTCD == "auclast") %>% .$PPORRES %>% 
            round(digits = 1)
    })
    
    output$pkauc1 <- renderText({
        # browser()
        reacauc()
    })
    
    
    dat_tumor <- reactive({
        
        # cyc <- input$cycle_tumor
        cyc <- 1
        
        amt1 <- input$amt_tumor1
        tau1 <- input$interval_tumor1
        n1 <- input$n_tumor1
        
        # amt2 <- input$amt_tumor2
        # tau2 <- input$interval_tumor2
        # n2 <- input$n_tumor2
        amt2 <- 0
        tau2 <- 0
        n2 <- 1
        
        cl <-input$cl_tumor
        v <- input$v_tumor
        ka <- input$ka_tumor
        fr <- input$f_tumor
        
        dos1 <- expand.ev(cmt = 1,
                          time = c(0,
                                   c(1:cyc) * (tau1*n1 + tau2*n2)),
                          amt = amt1,
                          ii = tau1,
                          addl = n1 - 1,
                          CL = cl,
                          V = v, 
                          KA = ka,
                          Fr = fr,
                          LAMDA0 = input$lamda0_tumor,
                          EMAX = input$emax_tumor,
                          EC50 = input$ec50_tumor,
                          # W0 = input$w0_tumor,
                          # FRAC = input$frac_tumor,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V = rnorm(1, 0, 0),
                          ETA_EMAX = rnorm(1, 0, 0),
                          ETA_EC50= rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos2 <- expand.ev(cmt = 1,
                          time = c(c(tau1*n1 + tau2*n2, 
                                     c(1:(cyc)) * (tau1*n1 + tau2*n2)) - tau2*n2)[-1],
                          amt = amt2,
                          ii = tau2,
                          addl = n2 - 1,
                          CL = cl,
                          V = v, 
                          KA = ka,
                          Fr = fr,
                          LAMDA0 = input$lamda0_tumor,
                          EMAX = input$emax_tumor,
                          EC50 = input$ec50_tumor,
                          # W0 = input$w0_tumor,
                          # FRAC = input$frac_tumor,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V = rnorm(1, 0, 0),
                          ETA_EMAX = rnorm(1, 0, 0),
                          ETA_EC50= rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos <- expand.ev(cmt = 1,
                         time = 0,
                         amt = 0,
                         ii = 24,
                         addl = 100,
                         CL = cl,
                         V = v, 
                         KA = ka,
                         Fr = fr,
                         LAMDA0 = input$lamda0_tumor,
                         EMAX = input$emax_tumor,
                         EC50 = input$ec50_tumor,
                         # W0 = input$w0_tumor,
                         # FRAC = input$frac_tumor,
                         ETA_CL = rnorm(1, 0, 0), 
                         ETA_V = rnorm(1, 0, 0),
                         ETA_EMAX = rnorm(1, 0, 0),
                         ETA_EC50= rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 2)
        
        union_all(dos1, dos2) %>% 
            union_all(dos) %>%
            arrange(ID, time, addl, ii)
        
    })
    
    output$distPlotTumorPK <- renderPlot({
        
        # cyc <- input$cycle_tumor
        cyc <- 1
        
        amt1 <- input$amt_tumor1
        tau1 <- input$interval_tumor1
        n1 <- input$n_tumor1
        
        # amt2 <- input$amt_tumor2
        # tau2 <- input$interval_tumor2
        # n2 <- input$n_tumor2
        amt2 <- 0
        tau2 <- 0
        n2 <- 1
        
        mod_tumor %>%
            data_set(dat_tumor()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>%
            plot(PKCONC~time)
    })
    
    output$distPlotTumorPD <- renderPlot({
        # cyc <- input$cycle_tumor
        cyc <- 1
        
        amt1 <- input$amt_tumor1
        tau1 <- input$interval_tumor1
        n1 <- input$n_tumor1
        
        # amt2 <- input$amt_tumor2
        # tau2 <- input$interval_tumor2
        # n2 <- input$n_tumor2
        amt2 <- 0
        tau2 <- 0
        n2 <- 1
        
        mod_tumor %>%
            data_set(dat_tumor()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>%
            plot(TUMOR~time)
    })
    
    output$tgi <- renderText({
        # browser()
        
        # cyc <- input$cycle_tumor
        cyc <- 1
        
        amt1 <- input$amt_tumor1
        tau1 <- input$interval_tumor1
        n1 <- input$n_tumor1
        
        # amt2 <- input$amt_tumor2
        # tau2 <- input$interval_tumor2
        # n2 <- input$n_tumor2
        amt2 <- 0
        tau2 <- 0
        n2 <- 1
        
        # browser()
        
        out <- mod_tumor %>%
            data_set(dat_tumor()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) 
        
        tmp <- out@data
        
        t0 <- tmp %>% filter(ID == 2) %>% .$TUMOR %>% max()
        t <- tmp %>% filter(ID == 1) %>% .$TUMOR %>% max()
        
        # browser()
        
        round(abs(t - t0)/t0 * 100, digits = 0)
        
        # reacauc()
    })
    
    
    
    dat_biomarker <- reactive({
        
        cyc <- input$cycle_biomarker
        
        amt1 <- input$amt_biomarker1
        tau1 <- input$interval_biomarker1
        n1 <- input$n_biomarker1
        
        amt2 <- input$amt_biomarker2
        tau2 <- input$interval_biomarker2
        n2 <- input$n_biomarker2
        
        cl <-input$cl_biomarker
        v <- input$v_biomarker
        ka <- input$ka_biomarker
        fr <- input$f_biomarker
        
        dos1 <- expand.ev(cmt = 1,
                          time = c(0,
                                   c(1:cyc) * (tau1*n1 + tau2*n2)),
                          amt = amt1,
                          ii = tau1,
                          addl = n1 - 1,
                          CL = cl,
                          V = v, 
                          KA = ka,
                          Fr = fr,
                          EMAX = input$emax_biomarker,
                          EC50 = input$ec50_biomarker,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V = rnorm(1, 0, 0),
                          ETA_EMAX = rnorm(1, 0, 0),
                          ETA_EC50= rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos2 <- expand.ev(cmt = 1,
                          time = c(c(tau1*n1 + tau2*n2, 
                                     c(1:(cyc)) * (tau1*n1 + tau2*n2)) - tau2*n2)[-1],
                          amt = amt2,
                          ii = tau2,
                          addl = n2 - 1,
                          CL = cl,
                          V = v, 
                          KA = ka,
                          Fr = fr,
                          EMAX = input$emax_biomarker,
                          EC50 = input$ec50_biomarker,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V = rnorm(1, 0, 0),
                          ETA_EMAX = rnorm(1, 0, 0),
                          ETA_EC50= rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos <- expand.ev(cmt = 1,
                         time = 0,
                         amt = 0,
                         ii = 24,
                         addl = 100,
                         CL = cl,
                         V = v, 
                         KA = ka,
                         Fr = fr,
                         EMAX = input$emax_biomarker,
                         EC50 = input$ec50_biomarker,
                         ETA_CL = rnorm(1, 0, 0), 
                         ETA_V = rnorm(1, 0, 0),
                         ETA_EMAX = rnorm(1, 0, 0),
                         ETA_EC50= rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 2)
        
        union_all(dos1, dos2) %>% 
            union_all(dos) %>%
            arrange(ID, time, addl, ii)
    })
    
    output$distPlotEfficacyPKPDtime <- renderPlot({
        cyc <- input$cycle_biomarker
        
        amt1 <- input$amt_biomarker1
        tau1 <- input$interval_biomarker1
        n1 <- input$n_biomarker1
        
        amt2 <- input$amt_biomarker2
        tau2 <- input$interval_biomarker2
        n2 <- input$n_biomarker2
        
        mod_biomarker %>%
            data_set(dat_biomarker()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>%
            plot(PKCONC + EFFECT~time)
    })
    
    output$distPlotEfficacyPKPD <- renderPlot({
        
        cyc <- input$cycle_biomarker
        
        amt1 <- input$amt_biomarker1
        tau1 <- input$interval_biomarker1
        n1 <- input$n_biomarker1
        
        amt2 <- input$amt_biomarker2
        tau2 <- input$interval_biomarker2
        n2 <- input$n_biomarker2
        
        mod_biomarker %>%
            data_set(dat_biomarker()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>%
            plot(EFFECT~PKCONC)
    })
    
    output$inhi <- renderText({
        
        # cyc <- input$cycle_biomarker
        cyc <- 1
        
        amt1 <- input$amt_biomarker1
        tau1 <- input$interval_biomarker1
        n1 <- input$n_biomarker1
        
        # amt2 <- input$amt_biomarker2
        # tau2 <- input$interval_biomarker2
        # n2 <- input$n_biomarker2
        amt2 <- 0
        tau2 <- 0
        n2 <- 1
        
        out <- mod_biomarker %>%
            data_set(dat_biomarker()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) 
        
        tmp <- out@data
        
        # browser()
        
        eff <- tmp %>% filter(ID == 1, 
                              time < cyc * (tau1*n1 + tau2*(n2-1)),
                              time > cyc * (tau1*(n1-1) + tau2*(n2-1))) %>% .$EFFECT %>% max()
        
        round(abs(eff) * 100, digits = 0)
    })
    
    
    
    
    dat_anc <- reactive({
        
        cyc <- input$cycle_anc
        
        amt1 <- input$amt_anc1
        tau1 <- input$interval_anc1
        n1 <- input$n_anc1
        
        amt2 <- input$amt_anc2
        tau2 <- input$interval_anc2
        n2 <- input$n_anc2
        
        cl <-input$cl_anc
        v <- input$v_anc
        ka <- input$ka_anc
        fr <- input$f_anc
        
        dos1 <- expand.ev(cmt = 1,
                          time = c(0,
                                   c(1:cyc) * (tau1*n1 + tau2*n2)),
                          amt = amt1,
                          ii = tau1,
                          addl = n1 - 1,
                          CL = cl,
                          V = v, 
                          KA = ka,
                          Fr = fr,
                          SLOPE = input$slope_anc,
                          # CIRC0 = input$circ0_anc,
                          GAMMA = input$gamma_anc,
                          MTT = input$mtt_anc,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos2 <- expand.ev(cmt = 1,
                          time = c(c(tau1*n1 + tau2*n2, 
                                     c(1:(cyc)) * (tau1*n1 + tau2*n2)) - tau2*n2)[-1],
                          amt = amt2,
                          ii = tau2,
                          addl = n2 - 1,
                          CL = cl,
                          V = v, 
                          KA = ka,
                          Fr = fr,
                          SLOPE = input$slope_anc,
                          # CIRC0 = input$circ0_anc,
                          GAMMA = input$gamma_anc,
                          MTT = input$mtt_anc,
                          ETA_CL = rnorm(1, 0, 0), 
                          ETA_V = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 1)
        
        dos <- expand.ev(cmt = 1,
                         time = 0,
                         amt = 0,
                         ii = 24,
                         addl = 100,
                         CL = cl,
                         V = v, 
                         KA = ka,
                         Fr = fr,
                         SLOPE = input$slope_anc,
                         # CIRC0 = 0.5,
                         GAMMA = input$gamma_anc,
                         MTT = input$mtt_anc,
                         ETA_CL = rnorm(1, 0, 0), 
                         ETA_V = rnorm(1, 0, 0)) %>% 
            mutate(dose = amt, ID = 2)
        
        
        union_all(dos1, dos2) %>% 
            union_all(dos) %>%
            arrange(ID, time, addl, ii)
        
    })
    
    output$distPlotSafetyPKPD <- renderPlot({
        cyc <- input$cycle_anc
        
        amt1 <- input$amt_anc1
        tau1 <- input$interval_anc1
        n1 <- input$n_anc1
        
        amt2 <- input$amt_anc2
        tau2 <- input$interval_anc2
        n2 <- input$n_anc2
        
        mod_anc %>%
            data_set(dat_anc()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) %>%
            plot(PKCONC + ANC~time)
    })
    
    output$nadir <- renderText({
        # browser()
        cyc <- input$cycle_anc
        
        amt1 <- input$amt_anc1
        tau1 <- input$interval_anc1
        n1 <- input$n_anc1
        
        amt2 <- input$amt_anc2
        tau2 <- input$interval_anc2
        n2 <- input$n_anc2
        
        # browser()
        
        out <- mod_anc %>%
            data_set(dat_anc()) %>%
            mrgsim(end = cyc * (tau1*n1 + tau2*n2),
                   delta = 0.05) 
        
        tmp <- out@data
        
        nadir <- tmp %>% filter(ID == 1) %>% .$ANC %>% min()
        
        # browser()
        
        round(nadir, digits = 2)
        
        # reacauc()
    })
    
    
})

