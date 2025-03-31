functions  {
    int[] find_start_end(int[] param_count, int t){

        int res[2];
        res[1] = 1;
        res[2] = 0;

      if (size(param_count) < t) {
        reject("selected layer index (for pair index) larger than the total number of layers; found layer_t (h) =", t);
      }
        //for starting index
        for (i in 1:(t-1)) {
          res[1] = res[1] + param_count[i];
        }

        //for ending index
        for (i in 1:t) {
          res[2] = res[2] + param_count[i];
        }
        return res;
      }
  }
data {
  int<lower=0> prior_sim; //boolean for whether to simulate from priors (1) or fit the actual model (0)
  int<lower=0> network_sim; //boolean for whether to draw dyadic outcomes from the posterior, 1 yes and 0 no
  int<lower=0> L; //number of groups/repeated obs of the multiplex network
  int<lower=0> n[L]; //number of actors for each of the L groups
  int<lower=0> N[L]; // number of total dyads for  [... truncated]
2: In file(fname, "rt") : expanded path length 12158 would be too long for
functions  {
    int[] find_start_end(int[] param_count, int t){

        int res[2];
        res[1] = 1;
        res[2] = 0;

      if (size(param_count) < t) {
        reject("selected layer index (for pair index) larger than the total number of layers; found layer_t (h) =", t);
      }
        //for starting index
        for (i in 1:(t-1)) {
          res[1] = res[1] + param_count[i];
        }

        //for ending index
        for (i in 1:t) {
          res[2] = res[2] + param_count[i];
        }
        return res;
      }
  }
data {
  int<lower=0> prior_sim; //boolean for whether to simulate from priors (1) or fit the actual model (0)
  int<lower=0> network_sim; //boolean for whether to draw dyadic outcomes from the posterior, 1 yes and 0 no
  int<lower=0> L; //number of groups/repeated obs of the multiplex network
  int<lower=0> n[L]; //number of actors for each of the L groups
  int<lower=0> N[L]; // number of total dyads for  [... truncated]
3: In file(fname, "rt") : expanded path length 12158 would be too long for
functions  {
    int[] find_start_end(int[] param_count, int t){

        int res[2];
        res[1] = 1;
        res[2] = 0;

      if (size(param_count) < t) {
        reject("selected layer index (for pair index) larger than the total number of layers; found layer_t (h) =", t);
      }
        //for starting index
        for (i in 1:(t-1)) {
          res[1] = res[1] + param_count[i];
        }

        //for ending index
        for (i in 1:t) {
          res[2] = res[2] + param_count[i];
        }
        return res;
      }
  }
data {
  int<lower=0> prior_sim; //boolean for whether to simulate from priors (1) or fit the actual model (0)
  int<lower=0> network_sim; //boolean for whether to draw dyadic outcomes from the posterior, 1 yes and 0 no
  int<lower=0> L; //number of groups/repeated obs of the multiplex network
  int<lower=0> n[L]; //number of actors for each of the L groups
  int<lower=0> N[L]; // number of total dyads for  [... truncated]
4: In file(fname, "rt") : cannot open file 'functions  {
    int[] find_start_end(int[] param_count, int t){

        int res[2];
        res[1] = 1;
        res[2] = 0;

      if (size(param_count) < t) {
        reject("selected layer index (for pair index) larger than the total number of layers; found layer_t (h) =", t);
      }
        //for starting index
        for (i in 1:(t-1)) {
          res[1] = res[1] + param_count[i];
        }

        //for ending index
        for (i in 1:t) {
          res[2] = res[2] + param_count[i];
        }
        return res;
      }
  }
data {
  int<lower=0> prior_sim; //boolean for whether to simulate from priors (1) or fit the actual model (0)
  int<lower=0> network_sim; //boolean for whether to draw dyadic outcomes from the posterior, 1 yes and 0 no
  int<lower=0> L; //number of groups/repeated obs of the multiplex network
  int<lower=0> n[L]; //number of actors for each of the L groups
  int<lower=0> N[L]; // number of total dyads for each of the L groups
  int<lowe [... truncated]
Error in get_model_strcode(file, model_code) : 
  cannot open model file "functions  {
    int[] find_start_end(int[] param_count, int t){

        int res[2];
        res[1] = 1;
        res[2] = 0;

      if (size(param_count) < t) {
        reject("selected layer index (for pair index) larger than the total number of layers; found layer_t (h) =", t);
      }
        //for starting index
        for (i in 1:(t-1)) {
          res[1] = res[1] + param_count[i];
        }

        //for ending index
        for (i in 1:t) {
          res[2] = res[2] + param_count[i];
        }
        return res;
      }
  }
data {
  int<lower=0> prior_sim; //boolean for whether to simulate from priors (1) or fit the actual model (0)
  int<lower=0> network_sim; //boolean for whether to draw dyadic outcomes from the posterior, 1 yes and 0 no
  int<lower=0> L; //number of groups/repeated obs of the multiplex network
  int<lower=0> n[L]; //number of actors for each of the L groups
  int<lower=0> N[L]; // number of total dyads for each of the L groups