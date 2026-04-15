#include <Rcpp.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include "pValueMethods.h"
//#include "malloc.h"

#include "Eigen/Dense"

using namespace Rcpp;
using Eigen::MatrixXd;

// Function annex to mcc, daudin, karlin, maxPartialSumD to :
//   *) check parameters consistency
//   *) allow to pass score values as vector (instead of just (min, max)) of the same size as score_probabilities
//   *) allow to pass score values as names(score_probabilities)
bool checkKDMparameters(int local_score, int sequence_length, NumericVector &score_probabilities, 
                        Rcpp::IntegerVector &sequence_min, Rcpp::IntegerVector &sequence_max,
                        Rcpp::IntegerVector &score_values) { 
  int scoremin ;
  int scoremax ;
  bool isnullseqmin = (sequence_min.size() == 0) ;
  bool isnullseqmax = (sequence_max.size() == 0) ;
  bool isnullscore = (score_values.size() == 0) ;
  bool isNoNameProb = (score_probabilities.attr("names")==R_NilValue) ;
  // bool isNoNameProb = (score_probabilities.attr("names").isUsable()) ;
  
  // if (isnullseqmin) sequence_min.push_front(R_PosInf) ;
  // if (isnullseqmax) sequence_max.push_front(R_NegInf) ;
  
  if(local_score<0)
    stop("[Invalid Input] local score must be positive.");
  if(sequence_length<= 0 )
    stop("[Invalid Input] sequence length must be positive.");
  double epsilon = 1e-12; // precision to compare the sum of the probability vector to 1.0
  if (fabs(sum(score_probabilities)-1.0)>epsilon)
    stop("[Invalid Input] score_probabilities must sum to 1.0.");

  if ((isnullseqmin && isnullseqmax && isnullscore && isNoNameProb) ||
      (!isnullseqmin && isnullseqmax) ||
      (isnullseqmin && !isnullseqmax)
  )
    stop("[Invalid input] Missing either sequence_max and sequence_min OR score_values vector OR integer names to score_probabilities.") ;
  
  // Convert names of score_probabilities into score values if not specify otherwise
  // if (isnullscore && isnullseqmin && !isNoNameProb) { 
  if (isnullscore && !isNoNameProb) { 
      CharacterVector ch = score_probabilities.names();
      if (ch.size() != score_probabilities.size()) {
        stop("[ERROR : Invalid Input] Either m is a matrix  with score values as rownames or specify score values via score_values parameter");
      } else { 
        // Initialisation of score values
        IntegerVector out(ch.size());
        std::transform(ch.begin(), ch.end(), out.begin(), std::atoi);
        score_values = out;
        // XXX TODO gestion erreur ou les rownames ne sont pas numeric (dans ce cas atoi renvoie '0')
      }
    }

  scoremin = isnullseqmin?min(score_values):sequence_min[0] ;
  scoremax = isnullseqmax?max(score_values):sequence_max[0] ;

  if (!isnullscore) {
    if (!isnullseqmin && (sequence_min[0] != min(score_values))) 
      stop("[Invalid input] sequence_min != min(score_values). Note : sequence_min and sequence_max are not required when score_values is given OR score_probabilities has scores as names.") ;
    if (!isnullseqmax && (sequence_max[0] != max(score_values))) 
      stop("[Invalid input] sequence_max != max(score_values). Note : sequence_min and sequence_max are not required when score_values is given OR score_probabilities has scores as names.") ;
    if (!isNoNameProb)
      warning("[Parameter specification] Both score_values and names(score_probabilities) are set. score_values is used.") ;
  }
  if (!isNoNameProb)
    if (!isnullseqmin)
      warning("[Parameter specification] Both sequence_min, sequence_max AND names(score_probabilities) are set. sequence_min, sequence_max are used.") ;
  
  // At this point, scoremin (resp. scoremax) contains the min (resp. max) of the scores in all cases
  
  if (scoremax<=0)
    stop("[Invalid Input] sequence_max must be positive.");
  if (scoremin>=0)
    stop("[Invalid Input] sequence_min must be negative.");
  if (!isnullscore && (score_probabilities.size()!= score_values.size()))
    stop("[Invalid Input] score probability distribution must contain as much elements as score_values.");
  if (!(isnullseqmin && isnullseqmax) && (score_probabilities.size() != scoremax - scoremin + 1))
    stop("[Invalid Input] score probability distribution must contain as much elements as the range from sequence_min to sequence_max. Note : instead of sequence_min and sequence_max, you can also use score_values vector, OR use names(score_probabilities) as score values.");

    // Remove element with zero probability at the beginning then at the end.
  unsigned int i = 0;
  while(std::fabs(score_probabilities[i])<1e-16) {
    scoremin++;
    score_probabilities = tail(score_probabilities,score_probabilities.length()-1);
    if(!isnullscore || !isNoNameProb)
      score_values = tail(score_values, score_values.length()-1) ;
  }
  i = score_probabilities.length()-1;
  while(std::fabs(score_probabilities[i])<1e-16) {
    scoremax--;
    score_probabilities = head(score_probabilities,score_probabilities.length()-1);
    if(!isnullscore || !isNoNameProb)
      score_values = head(score_values, score_values.length()-1) ;
    i--;
  }
  // Update sequence_min sequence_max
  if (isnullseqmin)
    sequence_min.push_front(scoremin) ;
  else
    sequence_min[0] = scoremin ;
  if (isnullseqmax)
    sequence_max.push_front(scoremax) ;
  else
    sequence_max[0] = scoremax ;

    // In case of score_values given, fill score_probabilities with 0. where no score is given
  if((!isnullscore || !isNoNameProb) &&  (score_values.size() != scoremax - scoremin + 1)) {  // otherwise nothing to do because no "holes" in score values
    // 1. Initialise un vecteur de 0.0 de la taille nécessaire en comblant les trous
    NumericVector res_prob (scoremax - scoremin + 1, 0.0) ; 
    // 2. Calcul les index des positions des scores
    for (int i = 0 ; i<score_values.size() ; i++)
      res_prob[score_values[i] - scoremin] = score_probabilities[i] ;
    score_probabilities = res_prob; 
  }

  return(true) ;
}

//' @description Calculates the exact p-value in the identically and independently distributed of a given local score, a sequence length that 'must not be too large' and for a given score distribution
//' @title Daudin [p-value] [iid]
//' @return A double representing the probability of a local score as high as the one given as argument. 
//' @param local_score the observed local score
//' @param sequence_length length of the sequence
//' @param score_probabilities the probabilities for each score from lowest to greatest (Optionnaly with scores as names)
//' @param sequence_min minimum score (optional if \code{score_values} OR names(score_probabilities) is defined)
//' @param sequence_max maximum score (optional if \code{score_values} OR names(score_probabilities) is defined)
//' @param score_values vector of integer score values, associated to score_probabilities  (optional if
//' \code{sequence_min} and \code{sequence_max} OR names(score_probabilities) are defined)
//' @details 
//' Either \code{sequence_min} and \code{sequence_max} are specified as input, OR all possible score values in 
//' \code{score_values} vector ; one of this choice is required. <cr>
//' Small in this context depends heavily on your machine. 
//' On a 3,7GHZ machine this means for daudin(1000, 5000, c(0.2, 0.2, 0.2, 0.1, 0.2, 0.1), -2, 3)
//' an execution time of ~2 seconds. This is due to the calculation method using matrix exponentiation
//' which takes times. The size of the matrix of the exponentiation is equal to a+1 with a the
//' local score value. The matrix must be put at the power n, with n the sequence length.
//' Moreover, it is known that the local score value is expected to be in mean of order log(n).
//' @examples 
//' p1 <- daudin(local_score = 4, sequence_length = 50, 
//'        score_probabilities = c(0.2, 0.3, 0.1, 0.2, 0.1, 0.1), 
//'        sequence_min = -3, sequence_max = 2)
//' p2 <- daudin(local_score = 4, sequence_length = 50, 
//'        score_probabilities = c(0.2, 0.3, 0.1, 0.2, 0.1, 0.1), 
//'        score_values = as.integer(-3:2))
//' p1 == p2 # TRUE
//' 
//' prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
//' score_values <- which(prob != 0) - 6 # keep only non null probability scores
//' prob0 <- prob[prob != 0]             # and associated probability
//' p <- daudin(150, 10000, prob, sequence_min = -5, sequence_max =  5)
//' p0 <- daudin(150, 10000, prob0, score_values = score_values)
//' names(prob0) <- score_values
//' p1 <- daudin(150, 10000, prob0)
//' p == p0 # TRUE
//' p == p1 # TRUE
//' @export
// [[Rcpp::export]]
double daudin(int local_score, int sequence_length, NumericVector score_probabilities,
              Rcpp::Nullable<Rcpp::IntegerVector> sequence_min = R_NilValue, 
              Rcpp::Nullable<Rcpp::IntegerVector> sequence_max = R_NilValue,
              Rcpp::Nullable<Rcpp::IntegerVector> score_values = R_NilValue) {
  Rcpp::IntegerVector score_values_ ; //non nullable equivalent
  Rcpp::IntegerVector sequence_min_ ; //non nullable equivalent
  Rcpp::IntegerVector sequence_max_ ; //non nullable equivalent
  sequence_min_ = sequence_min.isUsable() ? IntegerVector(sequence_min) : IntegerVector::create();
  sequence_max_ = sequence_max.isUsable() ? IntegerVector(sequence_max) : IntegerVector::create();
  score_values_ = score_values.isUsable() ? IntegerVector(score_values) : IntegerVector::create();
  
  bool checkOK = checkKDMparameters(local_score, sequence_length, score_probabilities, sequence_min_, sequence_max_, score_values_) ;

  return calcul_daudin(local_score, sequence_length, as<std::vector<double>>(score_probabilities),
                       as<int>(sequence_min_), as<int>(sequence_max_));
}

//' @description \code{karlin} Calculates an approximated p-value of a given local score value and a long sequence length in the identically and independently distributed model for the sequence. See also \code{\link{mcc}} function for another approximated method in the i.i.d. model that improved the one given by \code{\link{karlin}} or \code{\link{daudin}} for exact calculation. \cr
//' \code{karlin_parameters} is a annex function returning the parameters \eqn{\lambda}, \eqn{K^+} and \eqn{K^*} defined in Karlin and Dembo (1992).
//' @details This method works the better the longer the sequence is. Important note : the calculus of the parameter of the distribution uses
//' the resolution of a polynome which is a function of the score distribution, of order max(score)-min(score). There exists only empirical methods to solve a polynome of order greater that 5
//' with no warranty of reliable solution.
//' The found roots are checked internally to the function and an error message is throw in case of inconsistent. In such case, you could try to change your score scheme (in case of discretization)
//' or use the function \code{\link{karlinMonteCarlo}} .
//' This function implements the formulae given in Karlin and Dembo (1992), page 115-6. As the score is discrete here (lattice score function),
//' there is no limit distribution of the local score with the size of the sequence, but an inferior and a superior bound
//' are given. The output of this function is conservative as it gives the upper bound for the p-value.
//' Notice the lower bound can easily be found as it is the same call of function with parameter value local_score+1.
//' @title Karlin [p-value] [iid]
//' @return A double representing the probability of a local score as high as the one given as argument
//' @inheritParams daudin
//' @seealso \code{\link{mcc}}, \code{\link{daudin}}, \code{\link{karlinMonteCarlo}}, \code{\link{monteCarlo}}
//' @examples 
//' karlin(150, 10000, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -5, 5)
//' p1 <- karlin(local_score = 15, sequence_length = 5000, 
//'        score_probabilities = c(0.2, 0.3, 0.1, 0.2, 0.1, 0.1), 
//'        sequence_min = -3, sequence_max = 2)
//' p2 <- karlin(local_score = 15, sequence_length = 5000, 
//'        score_probabilities = c(0.2, 0.3, 0.1, 0.2, 0.1, 0.1), 
//'        score_values = -3:2)
//' p1 == p2 # TRUE
//' 
//' prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
//' score_values <- which(prob != 0) - 6 # keep only non null probability scores
//' prob0 <- prob[prob != 0]             # and associated probability
//' p <- karlin(150, 10000, prob, sequence_min = -5, sequence_max =  5)
//' p0 <- karlin(150, 10000, prob0, score_values = score_values) 
//' names(prob0) <- score_values
//' p1 <- karlin(150, 10000, prob0)
//' p == p0 # TRUE
//' p == p1 # TRUE
//' @export
// [[Rcpp::export]]
double karlin(int local_score, int sequence_length, NumericVector score_probabilities,
              Rcpp::Nullable<Rcpp::IntegerVector> sequence_min = R_NilValue, 
              Rcpp::Nullable<Rcpp::IntegerVector> sequence_max = R_NilValue,
              Rcpp::Nullable<Rcpp::IntegerVector> score_values = R_NilValue) {
  Rcpp::IntegerVector score_values_ ; //non nullable equivalent
  Rcpp::IntegerVector sequence_min_ ; //non nullable equivalent
  Rcpp::IntegerVector sequence_max_ ; //non nullable equivalent
  sequence_min_ = sequence_min.isUsable() ? IntegerVector(sequence_min) : IntegerVector::create();
  sequence_max_ = sequence_max.isUsable() ? IntegerVector(sequence_max) : IntegerVector::create();
  score_values_ = score_values.isUsable() ? IntegerVector(score_values) :IntegerVector::create();
      
  bool checkOK = checkKDMparameters(local_score, sequence_length, score_probabilities, sequence_min_, sequence_max_, score_values_) ;

  double esp = 0.0 ;
  int min_score = as<int>(sequence_min_) ;
  int max_score = as<int>(sequence_max_) ;
  for (int i=min_score ; i<=max_score ; i++) {
    esp += i*score_probabilities[i-min_score] ;
  }
  if(esp>=0.0)
    stop("[Invalid Input] Score expectation must be strictly negative.");

    double p = calcul_karlin(local_score, as<std::vector<double>>(score_probabilities), as<int>(sequence_max_), -as<int>(sequence_min_), (long)sequence_length);
  if (fabs(p+1.0)<1e-10) //p==-1.0 in case of error in calcul_karlin (polynomial roots problem)
    stop("karlin() function cannot be used in your case due to numerical instability (polynomial roots solver). Check the documentation of 'karlin()' for details.\n You could try to change your scoring discretisation step or use karlinMonteCarlo()");
  if (fabs(p+2.0)<1e-10) //p==-2.0 in case of error in calcul_karlin 
    stop("ERROR karlin: u and/or v are not compatible with the size of 'distribution'");
  return p;
  }

// return a vector containing K_star, K+ and Lambda of Karlin distribution (voir p. 115-116 karlin et Dembo 1992)
//' @rdname karlin
//' @export
// [[Rcpp::export]]
NumericVector karlin_parameters(NumericVector score_probabilities, 
                                Rcpp::Nullable<Rcpp::IntegerVector> sequence_min = R_NilValue, 
                                Rcpp::Nullable<Rcpp::IntegerVector> sequence_max = R_NilValue,
                                Rcpp::Nullable<Rcpp::IntegerVector> score_values = R_NilValue) {
  Rcpp::IntegerVector score_values_ ; //non nullable equivalent
  Rcpp::IntegerVector sequence_min_ ; //non nullable equivalent
  Rcpp::IntegerVector sequence_max_ ; //non nullable equivalent
  sequence_min_ = sequence_min.isUsable() ? IntegerVector(sequence_min) : IntegerVector::create();
  sequence_max_ = sequence_max.isUsable() ? IntegerVector(sequence_max) : IntegerVector::create();
  score_values_ = score_values.isUsable() ? IntegerVector(score_values) : IntegerVector::create();
  
  bool checkOK = checkKDMparameters(1, 1, score_probabilities, sequence_min_, sequence_max_, score_values_) ;

  double esp = 0.0 ;
  int min_score = as<int>(sequence_min_) ;
  int max_score = as<int>(sequence_max_) ;
  for (int i=min_score ; i<=max_score ; i++) {
    esp += i*score_probabilities[i-min_score] ;
  }
  if(esp>=0.0)
    stop("[Invalid Input] Score expectation must be strictly negative.");
  
    
  NumericVector res = Rcpp::wrap(calcul_karlin_parameters(as<std::vector<double>>(score_probabilities), max_score, -min_score)) ;
  res.names() = CharacterVector({"K_star", "K_plus", "lambda"}) ;
  double epsilon = 1e-10 ;
  if (std::fabs(res[0] + 1.0)<epsilon) // ERROR in polynomial roots finding
    stop("karlin_parameters() function cannot be used in your case due to numerical instability (polynomial roots solver). Check the documentation of 'karlin()' for details.\n You could try to change your scoring discretisation step or use karlinMonteCarlo()");

  return(res) ;
}

//' @description Calculates an approximated p-value for a given local score value and a medium to long sequence length in the identically and independently distributed model
//' @details This methods is actually an improved method of Karlin and produces more precise results. It should be privileged whenever possible. \cr
//' As with karlin, the method works the better the longer the sequence. Important note : the calculus of the parameter of the distribution uses
//' the resolution of a polynome which is a function of the score distribution, of order max(score)-min(score). There exists only empirical methods to solve a polynome of order greater that 5
//' with no warranty of reliable solution.
//' The found roots are checked internally to the function and an error message is throw in case of inconsistency. In such case, you could try to change your score scheme (in case of discretization)
//' or use the function \code{\link{karlinMonteCarlo}} .
//' @title MCC [p-value] [iid]
//' @return A double representing the probability of a local score as high as the one given as argument
//' @inheritParams daudin
//' @seealso \code{\link{karlin}}, \code{\link{daudin}}, \code{\link{karlinMonteCarlo}}, \code{\link{monteCarlo}}
//' @examples 
//' mcc(40, 100, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -6, 4)
//' mcc(40, 10000, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -6, 4)
//' mcc(150, 10000, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -5, 5)
//' p1 <- mcc(local_score = 15, sequence_length = 5000, 
//'        score_probabilities = c(0.2, 0.3, 0.1, 0.2, 0.1, 0.1), 
//'        sequence_min = -3, sequence_max = 2)
//' p2 <- mcc(local_score = 15, sequence_length = 5000, 
//'        score_probabilities = c(0.2, 0.3, 0.1, 0.2, 0.1, 0.1), 
//'        score_values = -3:2)
//' p1 == p2 # TRUE
//' 
//' prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
//' score_values <- which(prob != 0) - 6 # keep only non null probability scores
//' prob0 <- prob[prob != 0]             # and associated probability
//' p <- mcc(150, 10000, prob, sequence_min = -5, sequence_max =  5)
//' p0 <- mcc(150, 10000, prob0, score_values = score_values)
//' names(prob0) <- score_values
//' p1 <- mcc(150, 10000, prob0)
//' p == p0 # TRUE
//' p == p1 # TRUE
//' @export
// [[Rcpp::export]]
double mcc(int local_score, int sequence_length, NumericVector score_probabilities, 
           Rcpp::Nullable<Rcpp::IntegerVector> sequence_min = R_NilValue, 
           Rcpp::Nullable<Rcpp::IntegerVector> sequence_max = R_NilValue,
           Rcpp::Nullable<Rcpp::IntegerVector> score_values = R_NilValue) {
  Rcpp::IntegerVector score_values_ ; //non nullable equivalent
  Rcpp::IntegerVector sequence_min_ ; //non nullable equivalent
  Rcpp::IntegerVector sequence_max_ ; //non nullable equivalent
  sequence_min_ = sequence_min.isUsable() ? IntegerVector(sequence_min) : IntegerVector::create();
  sequence_max_ = sequence_max.isUsable() ? IntegerVector(sequence_max) : IntegerVector::create();
  score_values_ = score_values.isUsable() ? IntegerVector(score_values) : IntegerVector::create();
  
  bool checkOK = checkKDMparameters(local_score, sequence_length, score_probabilities, sequence_min_, sequence_max_, score_values_) ;

  double esp = 0.0 ;
  int min_score = as<int>(sequence_min_) ;
  int max_score = as<int>(sequence_max_) ;
  for (int i=min_score ; i<=max_score ; i++) {
    esp += i*score_probabilities[i-min_score] ;
  }
  if(esp>=0.0)
    stop("[Invalid Input] Score expectation must be strictly negative.");
  
  double p = calcul_mcc(local_score, as<std::vector<double>>(score_probabilities), max_score, -min_score, (long)sequence_length);
  if (fabs(p+1.0)<1e-10) //p==-1.0 in case of error in calcul_mcc
    stop("mcc() function cannot be used in your case. Check the documentation of 'mcc()' for details.\n You could try to change your scoring discretisation step or use karlinMonteCarlo()");
  return p;
}

//' @description Calculates the distribution of the maximum of the partial sum process for a given value in the identically and independently distributed model
//' @details Implement the formula (4) of the article Mercier, S., Cellier, D., & Charlot, D. (2003). An improved approximation for assessing the statistical significance of molecular sequence features. Journal of Applied Probability, 40(2), 427-441. doi:10.1239/jap/1053003554 \cr
//' Important note : the calculus of the parameter of the distribution uses
//' the resolution of a polynome which is a function of the score distribution, of order max(score)-min(score). There exists only empirical methods to solve a polynome of order greater that 5
//' with no warranty of reliable solution.
//' The found roots are checked internally to the function and an error message is throw in case of inconsistency. 
//' @title Maximum of the partial sum [probability] [iid]
//' @return A double representing the probability of the maximum of the partial sum process equal to k
//' @param k value at which calculates the probability
//' @inheritParams daudin
//' @examples 
//' maxPartialSumd(10, c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02), -6, 4)
//' prob <- c(0.08, 0.32, 0.08, 0.00, 0.08, 0.00, 0.00, 0.08, 0.02, 0.32, 0.02)
//' score_values <- which(prob != 0) - 7 # keep only non null probability scores
//' prob0 <- prob[prob != 0]             # and associated probability
//' p <- maxPartialSumd(10, prob, sequence_min = -6, sequence_max =  4)
//' p0 <- maxPartialSumd(10,prob0, score_values = score_values)
//' p == p0 # TRUE
//' @export
// [[Rcpp::export]]
double maxPartialSumd(int k, NumericVector score_probabilities,
                      Rcpp::Nullable<Rcpp::IntegerVector> sequence_min = R_NilValue, 
                      Rcpp::Nullable<Rcpp::IntegerVector> sequence_max = R_NilValue,
                      Rcpp::Nullable<Rcpp::IntegerVector> score_values = R_NilValue) {
  Rcpp::IntegerVector score_values_ ; //non nullable equivalent
  Rcpp::IntegerVector sequence_min_ ; //non nullable equivalent
  Rcpp::IntegerVector sequence_max_ ; //non nullable equivalent
  sequence_min_ = sequence_min.isUsable() ? IntegerVector(sequence_min) : IntegerVector::create();
  sequence_max_ = sequence_max.isUsable() ? IntegerVector(sequence_max) : IntegerVector::create();
  score_values_ = score_values.isUsable() ? IntegerVector(score_values) : IntegerVector::create();
  
  bool checkOK = checkKDMparameters(1, 1, score_probabilities, sequence_min_, sequence_max_, score_values_) ;

  double esp = 0.0 ;
  int min_score = as<int>(sequence_min_) ;
  int max_score = as<int>(sequence_max_) ;
  for (int i=min_score ; i<=max_score ; i++) {
    esp += i*score_probabilities[i-min_score] ;
  }
  if(esp>=0.0)
    stop("[Invalid Input] Score expectation must be strictly negative.");
  if(k<0)
     stop("[Invalid Input] local score must be strictly positive.");

     double p = calcul_probMaxPartialSum(k, as<std::vector<double>>(score_probabilities), max_score, -min_score);
   if (fabs(p+1.0)<1e-10) //p==-1.0 in case of error in calcul_probMaxPartialSum
     stop("probMaxPartialSum() function cannot be used in your case. Check the documentation of 'probMaxPartialSum()' for details.\n You could try to change your scoring discretisation step or use karlinMonteCarlo()");
   return p;
 }


//' @description Calculates stationary distribution of markov transition matrix by use of eigenvectors of length 1
//' @title Stationary distribution [Markov chains]
//' @return A vector with the probabilities
//' @param m Transition Matrix [matrix object]
//' @examples 
//' B = t(matrix (c(0.2, 0.8, 0.4, 0.6), nrow = 2))
//' stationary_distribution(B)
//' @export
// [[Rcpp::export]]
NumericVector stationary_distribution(NumericMatrix m){
  // std::accumulate(x.begin(),x.end(), 0.0)
  float rowsum ;
  float compliance = 0.00001 ;
  int i, j ;
  for(i = 0; i< m.nrow(); i++){
    NumericMatrix::Row row_i = m(i,_) ;
    rowsum = 0.0 ;
    for ( j = 0; j<row_i.size();j++) {
      rowsum += row_i[j] ;
    }
    if ( fabs(rowsum-1.0) > compliance ) {
      String the_message = "[ERROR] Transition probability matrix is not stochastic (row sum not equal 1.). Sum of line " + std::to_string(i+1) + " equal " + std::to_string(rowsum) ;
      stop(the_message) ;
      //stop(std::strcat("[ERROR] Transition probability matrix is not stochastic (row sum not equal 1.)",to_string(rowsum))) ;
      //stop(strcat("[ERROR] Transition probability matrix is not stochastic (row sum not equal 1.)", rowsum);
    }
  }

  MatrixXd mm(m.nrow(), m.ncol());
  for(int i = 0; i< m.nrow(); i++){
    for(int j = 0; j< m.ncol(); j++)
      mm(i, j) = m(i, j);
  }
  
  // check that m is a stochastic matrix (row sum == 1 for each row)
  for(int i = 0; i< m.nrow(); i++){
    for(int j = 0; j< m.ncol(); j++)
      mm(i, j) = m(i, j);
  }
  std::vector<Eigen::VectorXcd> v = stationary_distribution(mm);
  NumericVector results(v[0].size());
  if(v.size()> 1){
    stop("Markov matrix is not irreductible (many eigenvalues == 1).");
    return results; //useless in case of 'stop' above
  }
  if(v.size() == 0){
    stop("no eigenvector found.");
    return results;  //useless in case of 'stop' above
  } else {
    for(int i = 0; i < v[0].size(); i++)
      results[i] = v[0](i).real();
    return results;
  }
}

// XXX Pour les problèmes de warning CRAN concernant la doc, voir l'astuce dans la deuxième réponse : https://stackoverflow.com/questions/26263479/function-param-default-value-stdvector-initialization-with-rcpp-and-c11
// Voir aussi : https://github.com/RcppCore/Rcpp/issues/190 et https://github.com/Rcpp11/attributes/issues/36
//' @description Calculates the exact p-value for short numerical Markov chains. Memory usage and time computation can be too large for a high local score value and high score range (see details).
//' @title Exact method for p-value [Markov chains]
//' @return A double representing the probability of a local score as high as the one given as argument
//' @param local_score Integer local score for which the p-value should be calculated
//' @param m Transition matrix [matrix object]. Optionality, rownames can be corresponding score values. m should be a transition matrix of an ergodic Markov chain.
//' @param sequence_length Length of the sequence
//' @param score_values A integer vector of sequence score values (optional). If not set, the rownames of m are used if they are numeric and set.
//' @param prob0 Vector of probability distribution of the first score of the sequence (optional). If not set, the stationnary distribution of m is used.
//' @details This method computation needs to allocate a square matrix of size local_score^(range(score_values)). This matrix is then exponentiated to sequence_length.
//' @seealso \code{\link{monteCarlo}}
//' @examples 
//' mTransition <- matrix(c(0.2, 0.3, 0.5, 0.3, 0.4, 0.3, 0.2, 0.4, 0.4), byrow = TRUE, ncol = 3)
//' scoreValues <- -1:1
//' initialProb <- stationary_distribution(mTransition)
//' exact_mc(local_score = 12, m = mTransition, sequence_length = 100, 
//'         score_values = scoreValues, prob0 = initialProb)
//' exact_mc(local_score = 150, m = mTransition, sequence_length = 1000, 
//'          score_values = scoreValues, prob0 = initialProb)
//' rownames(mTransition) <- scoreValues
//' exact_mc(local_score = 12, m = mTransition, sequence_length = 100, prob0 = initialProb)
//' # Minimal specification
//' exact_mc(local_score = 12, m = mTransition, sequence_length = 100)
//' # Score valeus with "holes"
//' scoreValues <- c(-2, -1, 2)
//' mTransition <- matrix(c(0.2, 0.3, 0.5, 0.3, 0.4, 0.3, 0.2, 0.4, 0.4), byrow = TRUE, ncol = 3)
//' initialProb <- stationary_distribution(mTransition)
//' exact_mc(local_score = 50, m = mTransition, sequence_length = 100, 
//'         score_values = scoreValues, prob0 = initialProb)
//' @export
// [[Rcpp::export]]
double exact_mc(int local_score, NumericMatrix m, int sequence_length,
                Nullable<NumericVector> score_values = R_NilValue, 
                Nullable<NumericVector> prob0 = R_NilValue){
//  double exact_mc(int local_score, NumericMatrix m, int sequence_length, NumericVector score_values = Rcpp::NumericVector::create(), NumericVector prob0 = Rcpp::NumericVector::create()){
    
  double eps = 1e-12 ;
  NumericVector score_values_ ;
  NumericVector prob0_ ;
  if(score_values.isUsable()) {
    score_values_=score_values;
  } else {
    score_values_ = NumericVector::create();
  }
  if(prob0.isUsable()) {
    prob0_=prob0;
  } else {
    prob0_ = NumericVector::create();;
  }
  if(local_score<0)
    stop("[Invalid Input] local score must be positive.");
  if(sequence_length<= 0 )
    stop("[Invalid Input] sequence length must be positive.");
   if (m.nrow()!= m.ncol())
     stop("[ERROR exact_mc : Invalid Input] m should be a square matrix");
   if ((score_values_.size() != 0) && (m.nrow()!= score_values_.size() || m.ncol()!= score_values_.size())) 
     stop("[ERROR exact_mc : Invalid Input] m should be a square matrix of size the length of score_values");
   if ((prob0_.size() != 0) && (m.nrow()!= prob0_.size() || m.ncol()!= prob0_.size())) 
     stop("[ERROR exact_mc : Invalid Input] prob0 size should be equal to the number of rows of m");
   if ((prob0_.size() != 0) && (fabs(sum(prob0_) - 1.0)> eps))
     stop("[ERROR exact_mc : Invalid Input] prob0 vector should sum to 1");
   for(NumericVector::iterator i = prob0_.begin(); i != prob0_.end(); ++i) {
     //if ((fabs(*i) > eps) || (fabs(*i-1.0) > eps)) {
       if ((*i < -eps) || (*i > 1.0+eps)) {
         stop("[ERROR exact_mc : Invalid Input] prob0 vector should contains values between 0 and 1");
     }
   }
   
   // Deal with default parameters of R call
   if (score_values_.size() == 0) { // utilisation de la valeur par défaut
     if ((m.attr("dimnames")==R_NilValue) || (rownames(m)==R_NilValue)){
       stop("[ERROR exact_mc : Invalid Input] Either m is a matrix  with score values as rownames or specify score values via score_values parameter");
     } else { 
       CharacterVector ch = rownames(m);
       //     if (Rf_isNull(ch)) {
       if (ch.size() != m.ncol()) {
         stop("[ERROR exact_mc : Invalid Input] Either m is a matrix  with score values as rownames or specify score values via score_values parameter");
       } else { 
         // Initialisation of score values
         NumericVector out(ch.size());
         std::transform(ch.begin(), ch.end(), out.begin(), std::atoi);
         score_values_ = out;
         // XXX TODO gestion erreur ou les rownames ne sont pas numeric
       }
     }
   }
   if(max(score_values_)<=0)
    stop("[Invalid Input] sequence_max must be positive.");
   if(min(score_values_)>=0)
    stop("[Invalid Input] sequence_min must be negative.");
   if(prob0_.size() == 0) { // utilisation de la valeur par défaut (donné par stationary_distribution(m))
     NumericVector stationnary = stationary_distribution(m);
     prob0_ = stationnary ;
   }
  MatrixXd mm(m.nrow(), m.ncol());
   Eigen::VectorXi scoreValues(m.nrow());
   Eigen::VectorXd probInitial(m.nrow());
  for(int i = 0; i< m.nrow(); i++){
     scoreValues(i) = score_values_(i) ;
     probInitial(i) = prob0_(i) ;
    for(int j = 0; j< m.ncol(); j++)
      mm(i, j) = m(i, j);
  }
   return mh_markov(local_score, mm, scoreValues, sequence_length, probInitial);
}
