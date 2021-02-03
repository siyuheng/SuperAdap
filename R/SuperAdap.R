#' One-sided test in a sensitivity analysis in matched observational studies
#' via the two-stage programming method
#'
#' \code{SuperAdap} is the main function for performing a one-sided test
#'    under the Rosenbaum bounds sensitivity analysis framework via the two-stage
#'    programming method, which is an adaptive approach with a nice asymptotic
#'    property called "super-adaptivity". For details of the two-stage programming
#'    method, see the paper "Increasing Power for Observational Studies of Aberrant
#'    Response: An Adaptive Approach" by Heng, Kang, Small, and Fogarty.
#'
#' @param Q A N times 2 matrix of scores of the two component sum test statistics,
#'          where N is the total number of units in the study. That is, the n-th
#'          entry of the first column of \code{Q} is the score of the first component test
#'          statistic of unit n and the n-th entry of the second column of \code{Q} is the
#'          score of the second component test statistic of unit n.
#' @param Z A N-length vector of the treatment indicators of all N units in the
#'          study: 1 if treated and 0 if not. The n-th entry of \code{Z} is the treatment
#'          indicator of unit n. Note that in each matched set there is one and only
#'          one treated unit.
#' @param index A N-length vector of the matched set indexes of all N units in
#'              the study, taking value from 1 to I, where I is the total number of
#'              matched sets. That is, the n-th entry of \code{index} is the matched
#'              set index of unit n.
#' @param alpha The level of the one-sided test.
#' @param Gamma The sensitivity parameter in the Rosenbaum bounds sensitivity analysis, which
#'              is a prespecified number that is greater than or equal to 1.
#' @param alternative The direction of the alternative in a one-sided test, can be either
#'                  "greater", i.e., greater than, or "less", i.e., less than.
#' @return An indicator of the null hypothesis to be rejected or not: "reject" or "failed to reject".
#' @examples
#' #We randomly generate a dataset with I matched sets along with
#' #the scores of the two component test statistics
#' I=200
#' n<-rep(0, I)
#' for (i in 1:I){
#'   n[i]=sample(c(2:6), size = 1)
#' }
#' N=sum(n)
#' index_1<-rep(0, N)
#' Z_1<-rep(0, N)
#' Q_1<-matrix(0, nrow = N, ncol = 2)
#' S<-rep(0, I)
#' for (i in 1:I){
#'   S[i]=sum(n[1:i])
#' }
#' index_1[1:S[1]]=1
#' Z_1[1]=1
#' Z_1[1:S[1]]=sample(Z_1[1:S[1]], size = S[1])
#' if (I>1){
#'   for (i in 2:I){
#'     index_1[(S[i-1]+1):S[i]]=i
#'     Z_1[(S[i-1]+1)]=1
#'     Z_1[(S[i-1]+1):S[i]]=sample(Z_1[(S[i-1]+1):S[i]], size = n[i])
#'   }
#' }
#' for (s in 1:N){
#'   if (Z_1[s]==1){
#'    Q_1[s, 1]=rnorm(1, mean = 0.5, sd=1)
#'    Q_1[s, 2]=rnorm(1, mean = 0.6, sd=1)
#'  }
#'  else {
#'    Q_1[s, 1]=rnorm(1, mean = 0.1, sd=1)
#'    Q_1[s, 2]=rnorm(1, mean = 0, sd=1)
#'  }
#' }
#' #We then run the SuperAdap function with the generated dataset
#' result_1=SuperAdap(Q = Q_1, Z = Z_1, index = index_1,
#'                       alpha = 0.05, Gamma = 1.5, alternative = "greater")
#' result_2=SuperAdap(Q = Q_1, Z = Z_1, index = index_1,
#'                       alpha = 0.05, Gamma = 5, alternative = "greater")
#' @importFrom stats optim qchisq qnorm runif sd uniroot
#' @export



SuperAdap<-function(Q, Z, index, alpha, Gamma, alternative){
  if (!requireNamespace("gurobi", quietly = TRUE)) {
    stop("Package \"gurobi\" is required for this function to work.")
  }
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop("Package \"mvtnorm\" is required for this function to work.")
  }
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package \"Matrix\" is required for this function to work.")
  }

  if (alternative=='less'){
    Q=-Q
  }

  finderFS = function(q, Co, alpha)
  {
    1 - mvtnorm::pmvnorm(upper = rep(q, 2), corr = matrix(c(1,Co, Co, 1), nrow = 2, ncol = 2))[1] - alpha
  }

  minCorFullMatching = function(Q, index, Gamma)
  {
    epsilon = 1e-5
    K = 2
    nostratum = length(unique(index))
    noIndiv = length(index) # Number of individuals

    CovPair_FullMatch = function(expyu, Q, index)
    {
      covariance = 0
      for(i in 1:nostratum)
      {
        ind = which(index == i)
        qIthStrat = Q[ind, ]
        expyuIthStrat = expyu[ind] #the members of expyu that correspond to the individuals in the i^th stratum
        crossTerm = sum(qIthStrat[,1]*qIthStrat[,2]*expyuIthStrat) / sum(expyuIthStrat)
        productOfExpectations = sum(qIthStrat[,1]*expyuIthStrat) * sum(qIthStrat[,2]*expyuIthStrat) / (sum(expyuIthStrat)^2)
        covariance = covariance + (crossTerm - productOfExpectations)
      }
      return (covariance)
    }

    #computes the pairwise correlation
    CorPairMin_FullMatch = function(expyu, Q, index)
    {
      sd1 = sqrt(CovPair_FullMatch(expyu, cbind(Q[,1], Q[,1]), index))
      sd2 = sqrt(CovPair_FullMatch(expyu, cbind(Q[,2], Q[,2]), index))
      corr.ret = CovPair_FullMatch(expyu, Q, index) / (sd1 * sd2)
      corr.ret
    }

    #outputs the gradient of the covariance function with respect to the expyui terms (updated)
    GradCovPairMin_FullMatch = function(expyu, Q, index)
    {
      grad = rep(0, noIndiv)
      for(i in 1:nostratum)
      {
        ind = which(index == i)
        qIthStrat = Q[ind, ]

        expyuIthStrat = expyu[ind]

        crossTerm = (qIthStrat[,1]*qIthStrat[,2]*sum(expyuIthStrat) - sum(qIthStrat[,1]*qIthStrat[,2]*expyuIthStrat)) / (sum(expyuIthStrat)^2)

        productTerm = ((qIthStrat[,1]*sum(qIthStrat[,2]*expyuIthStrat) + qIthStrat[,2]*sum(qIthStrat[,1]*expyuIthStrat)) * (sum(expyuIthStrat)^2) -
                         (sum(qIthStrat[,1]*expyuIthStrat)*sum(qIthStrat[,2]*expyuIthStrat)*2*sum(expyuIthStrat))) / (sum(expyuIthStrat)^4)

        grad[ind] = crossTerm - productTerm
      }
      return (grad)
    }

    #outputs the gradient of the pairwise correlation function with respect to the elements of expyu
    GradCorPairMin_FullMatch = function(expyu, Q, index)
    {
      sd1 = sqrt(CovPair_FullMatch(expyu, cbind(Q[,1], Q[,1]), index))
      sd2 = sqrt(CovPair_FullMatch(expyu, cbind(Q[,2], Q[,2]), index))

      productDerivative = (sd2 * GradCovPairMin_FullMatch(expyu, cbind(Q[,1], Q[,1]), index) / (2 * sd1)) + (sd1 * GradCovPairMin_FullMatch(expyu, cbind(Q[,2], Q[,2]), index) / (2 * sd2))
      numerator = (GradCovPairMin_FullMatch(expyu, Q, index) * sd1 * sd2) - (CovPair_FullMatch(expyu, Q, index) * productDerivative)
      grad = numerator / ((sd1 * sd2)^2)
      grad
    }


    #expyu = rep((1 + Gamma) / 2 ,noIndiv) #midpoint of feasible region
    expyu = runif(n = noIndiv, min = 1, max = Gamma) # random start point
    optim(par = expyu, fn = CorPairMin_FullMatch, gr = GradCorPairMin_FullMatch, Q=Q, index=index, method = "L-BFGS-B", lower = 1 + epsilon, upper = Gamma - epsilon, control = list(lmm = 5, pgtol = 1e-8, fnscale = 1e-4))$value
  }

  multipleComparisonsRoot = function(Gamma,index, Q, Z, alpha, alternative , critval)
  {
    ns = table(index)
    ms = table(index[Z==1])

    nostratum = length(unique(index))

    if(is.null(dim(Q)))
    {
      Q = t(t(Q))
    }
    K = ncol(Q)
    if(any(ms!=1))
    {
      stop("Strata must have either one treated and the rest controls")
    }
    if(any(Z!=0 & Z!=1))
    {
      stop("Treatment Vector (Z) Must be Binary")
    }
    for(i in 1:nostratum)
    {
      if(ms[i] > 1)
      {
        ind = which(index==i)
        Z[ind] = 1-Z[ind]
        for(k in 1:K)
        {
          qsum = sum(Q[ind,k])

          Q[ind,k] = qsum - Q[ind,k]
        }
      }
    }

    treatment = (Z==1)
    PObefore = Q
    K = ncol(PObefore)

    sort.new = function(x)
    {
      temp = sort(unique(x))
      new = 1:length(temp)
      ret = rep(0,length(x))
      for(i in new)
      {
        ret[x == temp[i]] = i
      }
      ret
    }

    if(length(alternative)==1)
    {
      alternative = rep(alternative, K)
    }

    if(any(alternative != "two.sided" & alternative != "greater" & alternative != "less" & alternative != "TS" & alternative != "G"& alternative != "L"))
    {
      stop("Alternative options are two.sided/TS, greater/G, or less/L")
    }
    alternative[alternative=="less"] = "L"
    alternative[alternative=="greater"] = "G"
    alternative[alternative=="two.sided"] = "TS"

    indTS = which(alternative == "TS")
    indG = which(alternative == "G")
    indL = which(alternative == "L")
    nTS = length(indTS)
    nG = length(indG)
    nL = length(indL)
    orderq = c(indTS, indG, indL)
    PO= PObefore[,orderq, drop = F]


    Gamma.vec = Gamma
    Reject = rep(0, length(Gamma.vec))



    sds = apply(PO, 2, sd)
    SDS = matrix(sds, nrow(PO), ncol(PO), byrow = T)
    Qmat= (PO/SDS)

    ns = table(index)
    ns.types = (ns-1)
    N.total = nrow(Qmat)
    NS = matrix(ns[index], N.total, K)
    nullexpec = colSums(Qmat/NS)


    treatment = (1*treatment==1)
    Tobs = apply(Qmat[treatment,,drop = F], 2, sum)
    Tobsmat = matrix(Tobs, N.total, K, byrow = T)

    #bigM = 5*max((Tobs-nullexpec)^2)
    Q = Qmat
    Q2 = Qmat^2
    if(is.null(critval))
    {
      kappa = c(rep(qchisq(1-alpha/K,1),nTS), rep(qchisq(1-2*alpha/K,1),nG+nL))
    }else{kappa = rep(critval, K)}
    kappaMat =  matrix(kappa, nrow(PO), ncol(PO), byrow = T)

    Plin = -2*Tobsmat*Q - kappaMat*Q2
    Plinoneside = -kappaMat*Q2
    U = matrix(0, N.total, length(Gamma))
    RHO = U
    row.ind = rep(0, N.total + (K)*(2*N.total+nostratum+1)+ 4*N.total)
    col.ind = row.ind
    values = row.ind
    b = rep(0, K*(nostratum+1)+nostratum+2*N.total)
    nvariables = K*(nostratum+1)+N.total+nostratum+(K-nTS)+1
    for(i in 1:nostratum)
    {
      ind = which(index==i)
      row.ind[ind] = rep(i, length(ind))
      col.ind[ind] = ind
      values[ind] = 1
      b[i] = 1
    }

    quadcon = vector("list", K)

    for(k in 1:K)
    {
      if(k <= nTS)
      {
        for(i in 1:nostratum)
        {
          ind = which(index==i)
          row.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = rep((k-1)*(nostratum+1) + nostratum+i, length(ind)+1)
          col.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = c(ind, N.total+i + (k-1)*(nostratum+1))
          values[(k-1)*(2*N.total+nostratum+1)+c(N.total + ind, 2*N.total + i)] = c(-sqrt(kappa[k])*Q[ind,k], 1)
        }


        row.ind[(k-1)*(2*N.total+nostratum+1)+(2*N.total+nostratum+1):(3*N.total+nostratum+1)] = (k-1)*(nostratum+1)+c(rep(2*nostratum+1, N.total+1))
        col.ind[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(1:N.total, (k-1)*(nostratum+1)+N.total+nostratum+1)
        values[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(-Q[,k], 1)
        b[((k-1)*(nostratum+1) + nostratum+1):((k-1)*(nostratum+1) + 2*nostratum+1)] = c(rep(0, nostratum+1))
        rowq = N.total + ((k-1)*(nostratum+1)+1):(k*(nostratum+1))
        colq = rowq
        valq = rep(1, length(rowq))
        quadcon[[k]] = list()
        quadcon[[k]]$Qc = Matrix::sparseMatrix(rowq, colq, x = valq, dims = c(nvariables, nvariables))
        qq= rep(0, nvariables)
        qq[1:N.total] = Plin[,k]
        qq[N.total+K*(nostratum+1)+1] = -1

        quadcon[[k]]$q = qq
        quadcon[[k]]$rhs = -Tobs[k]^2

      }

      if(k>nTS)
      {
        for(i in 1:nostratum)
        {
          ind = which(index==i)
          row.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = rep((k-1)*(nostratum+1) + nostratum+i, length(ind)+1)
          col.ind[(k-1)*(2*N.total+nostratum+1)+ c(N.total + ind, 2*N.total + i)] = c(ind, N.total+i + (k-1)*(nostratum+1))
          values[(k-1)*(2*N.total+nostratum+1)+c(N.total + ind, 2*N.total + i)] = c(-sqrt(kappa[k])*Q[ind,k], 1)
        }

        row.ind[(k-1)*(2*N.total+nostratum+1)+(2*N.total+nostratum+1):(3*N.total+nostratum+1)] = (k-1)*(nostratum+1)+c(rep(2*nostratum+1, N.total+1))
        col.ind[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(1:N.total, (k-1)*(nostratum+1)+N.total+nostratum+1)
        values[(k-1)*(2*N.total+nostratum+1)+ (2*N.total+nostratum+1):(3*N.total+nostratum+1)] = c(Q[,k], 1)
        b[((k-1)*(nostratum+1) + nostratum+1):((k-1)*(nostratum+1) + 2*nostratum+1)] = c(rep(0, nostratum), Tobs[k])
        rowq = c(N.total + ((k-1)*(nostratum+1)+1):(k*(nostratum+1)-1),N.total + K*(nostratum+1)+1+(k-nTS))
        colq = rowq
        valq = rep(1, length(rowq))
        quadcon[[k]] = list()
        quadcon[[k]]$Qc = Matrix::sparseMatrix(rowq, colq, x = valq, dims =  c(nvariables, nvariables))
        qq= rep(0, nvariables)
        qq[1:N.total] = Plinoneside[,k]
        qq[N.total+K*(nostratum+1)+1] = -1

        quadcon[[k]]$q = qq
        quadcon[[k]]$rhs = 0

      }



    }

    mm = N.total+ K*(2*N.total+nostratum+1)+1
    rr = K*(nostratum+1)+nostratum+1
    cc = N.total + K*(nostratum+1)+2
    bind = rr
    if(nG!=0)
    {
      genmax  = vector("list", nG)

      for(k in (nTS+1):(nTS+nG))
      {
        ktemp = k-nTS
        genmax[[ktemp]]$resvar = cc
        genmax[[ktemp]]$vars = (k-1)*(nostratum+1)+N.total+nostratum+1
        genmax[[ktemp]]$con = 0
        cc=cc+1
      }
    }
    if(nL!=0)
    {
      genmin = vector("list", nL)
      for(k in (nTS+nG+1):(nTS+nG+nL))
      {
        ktemp = k-nTS-nG
        genmin[[ktemp]]$resvar = cc
        genmin[[ktemp]]$vars = (k-1)*(nostratum+1)+N.total+nostratum+1
        genmin[[ktemp]]$con = 0
        cc=cc+1
      }
    }


    mmnext = mm
    ccnext = cc
    rrnext = rr


    for(ee in 1:length(Gamma.vec))
    {
      Gamma.sens = Gamma.vec[ee]
      mm = mmnext
      cc = ccnext
      rr = rrnext
      for(i in 1:nostratum)
      {
        ind = which(index == i)
        for(j in ind)
        {
          row.ind[c(mm, mm+1)] = rep(rr, 2)
          col.ind[c(mm, mm+1)] = c(j, cc)
          values[c(mm, mm+1)] = 	c(1, -Gamma.sens)
          row.ind[c(  mm+2, mm+3)] = rep(rr+1, 2)
          col.ind[c(mm+2, mm+3)] = c(j, cc)
          values[c(mm+2, mm+3)]= c(-1, 1)
          rr = rr+2
          mm = mm+4
        }
        cc = cc+1
      }


      const.dir = c(rep("=", length(b)-2*N.total), rep("<=", 2*N.total))
      model = list()
      model$A = Matrix::sparseMatrix(row.ind, col.ind, x=values)
      model$sense = const.dir
      model$quadcon = quadcon
      model$rhs = b
      model$lb = c(rep(0, N.total), rep(-Inf, K*(1+nostratum)), -1e-4, rep(-Inf, nL+nG), rep(0, nostratum))
      model$obj = c(rep(0, length(model$lb) - nostratum-1-nG-nL), 1, rep(0, nostratum+nG+nL))
      model$ub = c(rep(Inf, N.total), rep(Inf, K*(1+nostratum) + 1), rep(Inf, nL+nG), rep(Inf, nostratum))
      model$vtypes = c(rep("C", N.total+K*(nostratum+1)+1), rep("C", nL+nG), rep("C", nostratum))
      model$modelsense = "min"
      model$objcons = 0

      if(nG>0)
      {
        model$genconmax = genmax
      }

      if(nL > 0)
      {
        model$genconmin = genmin
      }
      solm = gurobi::gurobi(model, params = list(OutputFlag = 0, Cutoff = 1e-4))
      if (solm$status == "INF_OR_UNBD")
      {
        Reject[ee] = -1e-4 # all feasible solutions lead us to fail to reject the null hypothesis

      }else if(solm$status == "CUTOFF"){
        Reject[ee] = 1 + 1/(1+Gamma)

      }else{
        Reject[ee] = solm$objval
      }
    }
    Reject
  }

  Co=minCorFullMatching(Q, index, Gamma)
  cv = (uniroot(finderFS, c(0, qnorm(1-alpha/4)), Co, alpha)$root)
  res = multipleComparisonsRoot(Gamma, index, Q, Z, alternative = "G", critval = cv^2)
  if (res>=0){
    return("reject")
  }
  else {
    return("failed to reject")
  }
}


