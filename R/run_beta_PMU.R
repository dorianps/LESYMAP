#' run_beta_PMU 
#'
#' Run Permutations for SVR-LSM
#' This function permutes the labels and generates a random SVR model & beta weights
#' for each voxel location / iteration.The number of permutation is specified via nperm
#' 
#' @param lesmat matrix of voxels (columns) and subjects (rows).
#' @param behavior vector of behavioral scores.
#'
#' @return matrix with a random SVR model including beta weights
#'
#' @author Daniel Wiesen; based on SVR-LSM Scripts by Zhang et al., 2014
#'
#' @export


	run_beta_PMU <- function(lesmat, behavior, ori_beta_val, maskIdx, posIdx, negIdx, maskDim, betaScale) {

	#shuffle labels
			permBehavior = sample(behavior)
			behavior = permBehavior
	
	#generate model
	 svr = svm(x = lesmat, y = behavior, scale = FALSE, type = 'eps-regression', kernel = 'radial', gamma = 5, cost = 30, epsilon = 0.1, na.action = na.omit)
  
	#get weights
	w = t(svr$coefs) %*% svr$SV
    pmu_beta_map = as.vector(w*betaScale)
	
	    tmp_map = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
        tmp_map[maskIdx] = pmu_beta_map
        pmu_beta_map = as.vector(tmp_map[maskIdx])
		ori_beta_val = as.vector(ori_beta_val)
	    ori_beta_val = ori_beta_val[maskIdx]

		
	    map_count_tmp = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
	    map_count_pos_tmp = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))
	    map_count_neg_tmp = array(0, dim=c(maskDim[1],maskDim[2],maskDim[3]))

        map_count_tmp[maskIdx]= map_count_tmp[maskIdx] + 1*(abs(pmu_beta_map) > abs(ori_beta_val))
        map_count_pos_tmp[posIdx] = map_count_pos_tmp[posIdx] + 1*(pmu_beta_map[posIdx] > ori_beta_val[posIdx])
        map_count_neg_tmp[negIdx] = map_count_neg_tmp[negIdx] + 1*(pmu_beta_map[negIdx] < ori_beta_val[negIdx])
	 
	   output = list(map_count_tmp=map_count_tmp, map_count_pos_tmp=map_count_pos_tmp, map_count_neg_tmp=map_count_neg_tmp)

	 return(output) 

	 
	 }