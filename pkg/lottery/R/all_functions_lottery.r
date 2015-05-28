
#' Function modified from the untb package
#' code{untb2}
#' @param a,start Starting ecosystem; coerced to class census. Usually, pass an object of class count; see examples. To start with a monoculture of size 10, use start=rep(1,10) and to use start=1:10.
#' @param prob,prob.of.immigrate,prob.of.mutate Probability of “new” organism not being a descendent of an existing individual
#' @param D Number of organisms that die in each timestep
#' @param gens Number of generations to simulate
#' @param keep In function untb() Boolean with default FALSE meaning to return the system at the end of the simulation and TRUE meaning to return a matrix whose rows are the ecosystem at successive times
#' @param meta In function untb(), the metacommunity; coerced to a count object. Default of NULL means to use a “greedy” system in which every mutation gives rise to a new, previously unencountered species. This would correspond to an infinitely large, infinitely diverse, Hubbellian ecosystem (which is not too ridiculous an assumption for a small island near a large diverse mainland). In function select.immigrate(), a simplified representation of a metacommunity.
#' @param limit.sim Force de la limite à la similarité sous la forme d''une matrice de distance écologique entre individus
#' @param coeff.lim.sim Module l''écart entre espèce la plus pénalisée et la moins pénalisée
#' @param sigma Module le rapport entre compétition intra et interspécifique
#' @param hab.filt Traits sélectionnés sous la forme d''un vecteur correspondant au valeur individuel 
#' @param prob.mort Probabilité de mort des individus trop proches

@



untb2 <- function (start, prob = 0, D = 1, gens = 150, keep = FALSE, meta = NULL, 
                  limit.sim = NULL, coeff.lim.sim = 2, sigma=0.1, 
                  hab.filt = NULL, prob.mort = 0, method.dist="euclidean") {
#__________
#Create the metacommunity pool
  if (!is.null(meta)) {
    jj <- as.count(start)
    jj[] <- 0
    meta <- as.count(meta) + jj
    jj.meta <- list(spp = as.numeric(names(meta)), abundance = meta)
  }
  
  else {jj.meta <- NULL}
  
  if(!is.null(limit.sim)){
    if(!class(limit.sim)=="dist"){
      limit.sim <- dist(limit.sim, method=method.dist)
    }
  }
  
#__________
  if (is.census(start) | is.count(start)) {
      a <- as.integer(as.census(start))
  }  
  else {a <- start}

#__________
  n <- length(a)
  limit.sim.glob <- NULL;
  if(!is.null(meta) & !is.null(limit.sim))  {
    #__________
    # Limite à la similarité globale dans la métacommunauté
    limit.sim <- as.matrix(limit.sim)
    diag(limit.sim) <- NA
    limit.sim <- as.dist(limit.sim)
    limit.sim.glob <- NULL;
    #limit.sim.glob <- apply(limit.sim[jj.meta$spp,jj.meta$spp],2,function(x) sum((max(x[!is.na(x)])-x[!is.na(x)])*jj.meta$abundance[!is.na(x)])/sum(jj.meta$abundance[!is.na(x)]))
    #limit.sim.glob <- log(limit.sim.glob) 
  }
  
#__________  
  if (keep) {
    aa <- matrix(NA, gens, n)
    limit.sim.t <- c(); w.limit.sim.t <- c(); 
    for (i in 1:gens) {
      aa[i, ] <- a
      a <- select(a, D = D, prob = prob, meta = jj.meta,  prob.mort= prob.mort, 
                  limit.sim=limit.sim, coeff.lim.sim = coeff.lim.sim, 
                  limit.sim.glob=limit.sim.glob, sigma=sigma, hab.filt=hab.filt, 
                  first=first)
      if(!is.null(limit.sim)) {
        limit.sim.t <- c(limit.sim.t,a$limit.sim.t)
        w.limit.sim.t <- c(w.limit.sim.t,a$w.limit.sim.t)
      }
        a <- a$a
    } 
    return(list(aa=aa,limit.sim.t=limit.sim.t,w.limit.sim.t=w.limit.sim.t))
  }
    
  else {
    for (i in 1:gens) {
      if(i==1) {first = TRUE} else {first = FALSE}
      a <- select(a, D = D, prob = prob, meta = jj.meta, prob.mort= prob.mort, limit.sim=limit.sim, coeff.lim.sim = coeff.lim.sim, limit.sim.glob=limit.sim.glob, sigma=sigma, hab.filt=hab.filt, first=first)
      #if(!is.null(limit.sim)) print(a$limit.sim.t)
      a <- a$a
    }
  return(as.count(a))
  }
}



# Return "mutate" individuals if meta = NULL, else return immigrate individuals
select <- function (a, D = length(a), prob = 0, meta = NULL, prob.mort= prob.mort, 
                    limit.sim=NULL, limit.sim.glob=NULL, coeff.lim.sim = 2, 
                    sigma=0.1, hab.filt=NULL, first=TRUE) {
  if (is.null(meta)){
    return(select.mutate(a, D = D, prob.of.mutate = prob))
  } 
  else {
    return(select.immigrate(a, D = D, prob.of.immigrate = prob, meta = meta, 
                            prob.mort= prob.mort, limit.sim=limit.sim, 
                            limit.sim.glob=limit.sim.glob, 
                            coeff.lim.sim = coeff.lim.sim, sigma=sigma, 
                            hab.filt=hab.filt,first=first))
  }
}



select.mutate <- function (a, D = length(a), prob.of.mutate = 0) {
    n <- length(a)
    died <- sample(n, D, replace = TRUE)
    mutated <- runif(length(died)) < prob.of.mutate
    n1 <- sum(mutated)
    n2 <- sum(!mutated)
    a[died[!mutated]] <- sample(a, n2, replace = TRUE)
    a[died[mutated]] <- (1:n1) + max(a)
    return(list(a=a))
}




select.immigrate <- function (a, D = 1, prob.of.immigrate = 0, 
                              meta, prob.mort= 0, limit.sim.glob = NULL, 
                              limit.sim = NULL, coeff.lim.sim = 2, sigma = 0.1, 
                              hab.filt = NULL, first=TRUE){
  n <- length(a);
  a.init <- a;
  a.sp <- table(a);
#____________
#Mortality including potential limiting similarity
  if(is.null(limit.sim)) {
    died <- sample(n, D, replace = TRUE)
    a <- a[-died]; if(sum(is.na(a))!=0) {
      stop("Error: NA values in community composition (1)")
    }
  }
  #____________
  #Mortality effect of limiting similarity
  else{
    if(sum(!names(a.sp)%in%rownames(limit.sim))>0 | sum(!names(a.sp)%in%colnames(limit.sim))>0){
      stop("limit.sim: mismatch of species names")
    }
    if(length(a.sp)==1){
      a.sp = a.sp - D; #warning("One species remaining"); # Pour débuggage
      limit.sim.t <- a.sp; #limit.sim[names(a.sp)];
      w.limit.sim.t <- limit.sim.t;
      #print(names(a.sp))  # Pour débuggage
    }
    else{
      #limit.sim.t <- apply(limit.sim[names(a.sp),names(a.sp)],2,min)  # Distance au plus proche voisin
      # Option avec prise en compte des abondances
      # Mode de calcul analogue ? intensit? g?n?rale de la comp?tition interindividuelle, en n?ligeant la comp?tition intrasp?cifique
      #limit.sim.t <- apply(limit.sim[names(a.sp),names(a.sp)],2,function(x) sum((max(x[!is.na(x)])-x[!is.na(x)])*a.sp[!is.na(x)])/sum(a.sp[!is.na(x)]))
      #limit.sim.t <- log(limit.sim.t) # Pour une distribution ? peu pr?s normale
      # Mode de calcul analogue ? intensit? g?n?rale de la comp?tition interindividuelle, AVEC comp?tition intrasp?cifique
      limit.sim.t <- apply(limit.sim[names(a.sp),names(a.sp)],2, 
                           function(x) (a.sp[is.na(x)]+
                             sum(exp(-x[!is.na(x)]^2/(2*sigma^2))*a.sp[!is.na(x)])))  
      #plot(limit.sim.t,apply(limit.sim[names(a.sp),names(a.sp)],2,function(x) min(x,na.rm=T)))
      #died.sp <- rmultinom(1,D,a.sp*prob.mort*((coeff.lim.sim-1)*(limit.sim.t[names(a.sp)]-min(limit.sim.glob))/(max(limit.sim.glob)-min(limit.sim.glob))+1));
      #died.sp <- rmultinom(1,D,a.sp*prob.mort*((coeff.lim.sim-1)*limit.sim.t[names(a.sp)]+1));
      #died.sp <- rmultinom(1,D,a.sp*prob.mort*((coeff.lim.sim-1)*(max(limit.sim.t)-limit.sim.t)/(max(limit.sim.t)-min(limit.sim.t))+1));
      died.sp <- rmultinom(1, D, prob.mort*a.sp*((coeff.lim.sim-1)*limit.sim.t[names(a.sp)]+1));
     
      if(sum(is.na(died.sp))>0) {stop("Erreur: NA in died individuals")}
      #if(first) {print(coeff.lim.sim); first=F;}
      w.limit.sim.t <- weighted.mean(limit.sim.t,w=a.sp)
      limit.sim.t <- mean(limit.sim.t)
      if(any((as.matrix(a.sp)-as.matrix(died.sp))<0)) {
        stop("Mortality cannot be modeled under limiting similarity when D is too large")
      }
      a.sp <- as.matrix(a.sp)-as.matrix(died.sp);
    }
        a <- unlist(lapply(which(a.sp!=0),function(x) rep(x,a.sp[x])));
        a <- as.numeric(rownames(a.sp)[a]);  # Composition de la communaut? apr?s mortalit?
        if(sum(is.na(a))!=0) stop("Error: NA values in community composition (2)");
  }#end of mortality by limiting similarity   

  immigrated <- runif(D) < prob.of.immigrate; 
  n1 <- sum(immigrated);
  n2 <- sum(!immigrated);
  if (n1 > 0) {
#____________
#Immigration including potential habitat filtering    
    if(is.null(hab.filt)){
      a <- c(a, sample(meta$spp, n1, replace = TRUE, prob = meta$abundance))
      if(sum(is.na(a))!=0) stop("Error: NA values in immigrants");
    }
    else 
        {
          if(sum(!meta$spp%in%names(hab.filt))>0) stop("habitat filtering: mismatch of species names");
          a <- c(a, sample(meta$spp, n1, replace = TRUE, prob = meta$abundance*hab.filt[meta$spp]))
          if(sum(is.na(hab.filt[meta$spp]))!=0) stop("Error: NA values in habitat filtering")
        }
        if(sum(is.na(a))!=0) {print(a);stop("Error: NA values in community composition (3)");}
    }

  if (n2 > 0) {
    if(is.null(limit.sim)) {
      a <- c(a,sample(a.init, n2, replace = TRUE))
    }
    #a[died[!immigrated]] <- rmultinom(a.sp,n2,a.sp*(max(limit.sim)-limit.sim)/(sum(a.sp)*sum(max(limit.sim)-limit.sim)));} # Effet limiting similarity durant l'immigration    
    else {a <- c(a,sample(a.init, n2, replace = TRUE))} # Pas d'effet limiting similarity durant l'immigration
    if(sum(is.na(a))!=0) {
      print(n2); stop("Error: NA values in community composition (4)")
    }
  }
  if(!is.null(limit.sim)){
    return(list(a=a, limit.sim.t=limit.sim.t, w.limit.sim.t=w.limit.sim.t))
  }
  else {return(list(a=a))}
}





Adaptation de la fonction Matlab eponyme. \texttt{generneutral} génère une distribution de fréquence d'espèces pour une communauté locale neutre.\\

\textbf{Arguments:}

\begin{itemize}
  \item J:  de nbre d'individus dans la communauté
  \item theta : nbre fondamental de biodiversité
  \item m : taux d'immigration. Pour m=1, on a une méta-communauté
  \item T (option) : Pool des ancetres "forcé" par une communauté environnante T est une liste d'individus étiquetés par ancetres et especes (cf. res).
\end{itemize}

\textbf{Résultats:}

\begin{itemize}
  \item Esp : distribution de fréquence des espèces
  \item res : étiquetage des individus en ancetres et en espèces
  \item An : distribution de fréquence des ancetres
  \item Dp : fréquences ancetres x espèces
\end{itemize}


generneutral <- function(J, theta, m, Tind=NULL)  {
  A<-array(0,c(J,1))
  S<-array(0,c(J,1))

  if(m<1) {
    I<-m*(J-1)/(1-m)
    X<-runif(J)
  	R1<-I/(I+(1:J)-1)
  	mi1<-which(X<=R1)
  	mi2<-which(X>R1)
  } 
  else {mi1=1:J; mi2=c()}

  # mi1: ancetres extérieurs ; mi2: ancetres locaux
  AA<-length(mi1); # Nbre total d'ancetres extérieurs
  A[mi1]<-1:AA;

  if(!is.null(Tind))    # "forçage" par un pool d'ancetres imposé
  {
      Tind <- Tind[sample(nrow(Tind),J),];
      A[mi1] <-Tind[mi1,1];
      S[mi1] <- Tind[mi1,2];
  }
  
  else
  {
      Tnov <- novesp(J,theta,A[mi1]);
      A[mi1] <- Tnov[,1];
      S[mi1] <- Tnov[,2];
  }

  for(j in 1:length(mi2))
  {
      jj<-max(round(runif(1)*mi2[j]-1), 1);
      A[mi2[j]]<-A[jj];
      S[mi2[j]]<-S[jj];
  }

  res <- cbind(A,S);

  Dp <- table(A,S);
  An <- rowSums(Dp);
  Esp <- colSums(Dp);

  Esp=table(S);

  ret<-c();
  ret$Esp<-Esp;ret$res<-res
  ret$Dp<-Dp
  ret$An<-An;
  return(ret)
}


\texttt{novesp} est une fonction interne de generneutral. \texttt{novesp} identifie les "nouvelles" espèces (par mutations) dans une liste (B) de migrants. 
novesp <- function(J, theta, B)
{
  # Identifie les "nouvelles" espèces (mutations) dans une liste (B) de migrants
  # A noter la grande analogie avec la fct principale

  AA<-length(B);
  A2<-array(0,c(AA,1));
  S2<-array(0,c(AA,1));
  aa<-1; ss<-0;
  Y<-runif(AA);

  R2<-theta/(theta+(1:AA)-1);
  mia<-which(Y<=R2);
  mib<-which(Y>R2);

  A2[mia]<- B[mia];
  S2[mia]<- 1:length(mia);

  for(j in 1:length(mib))
  {
      jj<-max(round(runif(1)*mib[j]-1),1);
      A2[mib[j]]<-A2[jj];
      S2[mib[j]]<-S2[jj];
  }

  return(cbind(A2,S2))
}





\subsection{generneutral.multi}
generneutral.multi simule \texttt{n} communautés neutres partageant un pool de migrants.
  
\begin{itemize}
  \item J:  de nbre d'individus dans la communauté
  \item theta : nbre fondamental de biodiversité
  \item m : taux d'immigration
\end{itemize}


generneutral.multi <- function(n, J, theta, m, pool=NULL)
{
  if(length(m) == 1 & length(m) != n) {
    m <- rep(m, times= n)
  }
    
  # Création d'un ensemble de communautés neutres partageant un pool de migrants
  if(is.null(pool)){
    pool <- generneutral(J, theta, 1)$res
  }
  else if(is.null(ncol(pool))) {
    # Par défaut le code sp est dupliqué
    pool <- cbind(pool, pool)
  }
   
  tabl <- array(0, c(max(pool[,2]), n));
  for(i in 1:n) {
   comp <- generneutral(J, theta, m[i], Tind=pool)$Esp
   tabl[as.numeric(names(comp)), i] <- comp
  }
   
   return(tabl)
}





Rappel:
-  a, start : Starting ecosystem; coerced to class census. Usually, pass an object of class count; see examples. To start with a monoculture of size 10, use start=rep(1,10) and to use start=1:10.

- prob, prob.of.immigrate, prob.of.mutate  : Probability of “new” organism not being a descendent of an existing individual

- D  : Number of organisms that die in each timestep

- gens : Number of generations to simulate

- keep : In function untb() Boolean with default FALSE meaning to return the system at the end of the simulation and TRUE meaning to return a matrix whose rows are the ecosystem at successive times

- meta : In function untb(), the metacommunity; coerced to a count object. Default of NULL means to use a “greedy” system in which every mutation gives rise to a new, previously unencountered species. This would correspond to an infinitely large, infinitely diverse, Hubbellian ecosystem (which is not too ridiculous an assumption for a small island near a large diverse mainland). In function select.immigrate(), a simplified representation of a metacommunity.

- limit.sim :  force de la limite à la similarité sous la forme d'une matrice de distance écologique entre individus ou d'une matrice (un vecteur) de traits

- coeff.lim.sim : module l'écart entre espèce la plus pénalisée et la moins pénalisée

- sigma : module le rapport entre compétition intra et interspécifique

- hab.filt : traits sélectionnés sous la forme d'un vecteur correspondant au valeur individuel. Si il y a plusieurs habitats, hab.filt est une matrice dont chaque colonne correspond à un type d'habitat. 


- prob.mort : (probabilité de mort des individus trop proches) \\




\textbf{Autres argument de la fonction \texttt{SimForw}:}

- ncom: nombre de communauté à simuler

- val_I : liste des valeurs de I à tester

- nrep: Nombre de simulations par valeur de paramêtre 

- nb_hab: nombre d'habitat à simuler, doit être égale aux nombres de colonnes de hab.filt et doit être un diviseur de hab.filt.

- theta: nbre fondamental de biodiversité (par défault 50)
\\


\textbf{Résultats de la fonction \texttt{SimForw}:}

une liste contenant: 

-\$res : une matrice avec 6 colonnes indiquant le numéro du paramètre et sa valeur I, le numéro (le nom) de la communauté, l'espéce et son nb d'individus associé et enfin le nombre maximum d'espèces dans cette communauté.
    
-\$val_I : Le vecteur des valeurs de paramètres de I

-\$start : La communauté de départ 

-\$limit.sim : valeur de l'argument de limiting similarity

-\$hab.filt  : valeur de l'argument d'habitat filtering 



#Ci après la fonction \texttt{SimForw} simule des communautées neutres forward en faisant varier limmigration (paramètre I).

SimForw <- function(n= 10, nrep=5, val_I=c(1,10), start = NULL, prob = NULL, D = 1, 
                    gens = 200, keep = FALSE, meta = NULL, limit.sim = NULL, 
                    coeff.lim.sim=NULL, sigma = 0.1, hab.filt = NULL, 
                    prob.mort = 0.1, nb_hab=1, theta = 50) {
  
  is.wholenumber <- function(x, tol = .Machine$double.eps)  abs(x - round(x)) < tol
  if(!is.wholenumber(n/nb_hab)){stop("n/nb_hab need to be a whole number!")}
  if(!is.null(hab.filt)) {
    if(dim(as.matrix(hab.filt))[2]!=nb_hab){
      stop("the number of columns of hab.filt need to be equal to nb_hab")
    }
  }
  
  J <- length(start) / n #J est le nb total d'individus dans la région
  res.list <- list();res.list$res<-list()
  for(indimm in 1:length(val_I)){
    Imm <- val_I[indimm]*10;
    m <- Imm/(Imm+J-1);
    RES <- c(); effesp <- c(); PoolMetaN <- c()
    #**********************
    for(i in 1:nrep){
      meta <- generneutral(n*J, m=1,theta=theta)$Esp;
      MaxSp <- as.numeric(names(meta)[length(meta)]);
      Esp <- array(0, c(MaxSp, 1)); names(Esp) <- as.character(1:MaxSp);
      Esp[names(meta)] <- meta;
      PoolMetaN <- rbind(PoolMetaN, cbind(i*rep(1,MaxSp), Esp));
      #**********************
      for(h in 1:nb_hab){
        for(j in 1:(n/nb_hab)) {
          if(!is.null(hab.filt)){
            hab.filtering = as.vector(hab.filt[,h])
            names(hab.filtering) <- rownames(hab.filt)
          }
          loc <- untb(start=start, prob = m, D = D, gens = gens, keep = keep, meta = meta,
                      limit.sim = limit.sim, coeff.lim.sim=coeff.lim.sim, sigma = sigma,
                      hab.filt = hab.filtering, prob.mort = prob.mort) 
          if(keep) {
            loc <- loc$aa; loc<-loc[nrow(loc),];tab.sp <- table(loc)
          } 
          else {tab.sp <- loc}
          nesp = length(tab.sp);
          RES = rbind(RES, cbind(i*rep(1, nesp), val_I[indimm]*rep(1, nesp),
                                 j*rep(1, nesp), 
                                 as.numeric(rownames(tab.sp)), tab.sp,
                                 MaxSp*rep(1, nesp)));
          colnames(RES) <- c("Num_repet_param", "Value_I", "Num_Communauté", 
                             "Species", "Nb_ind_Sp" ,"Max_nb_Sp")
        }
      }
      print(paste(((1:length(val_I))[indimm])/length(val_I) * 100, "%", 
                  i/nrep * 100, "%", sep="  "))
    }
    res.list$res[[indimm]] <- RES
    res.list$val_I <- val_I
    res.list$start <- start
    res.list$limit.sim <- limit.sim
    res.list$hab.filt <- hab.filt
    res.list$call <- match.call()
  }
  return(res.list)
}
