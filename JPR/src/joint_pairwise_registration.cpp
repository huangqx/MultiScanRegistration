#include "joint_pairwise_registration.h"
#include "dynamic_linear_algebra.h"
#include "linear_algebra.h"
#include <algorithm>
#include <set>
using namespace std;

void ANNDataStructure::ClearANN() {
  if (dists_ != NULL) {
    delete[]dists_;
    dists_ = NULL;
  }
  if (nnIdx_ != NULL) {
    delete[]nnIdx_;
    nnIdx_ = NULL;
  }
  if (queryPt_ != NULL) {
    annDeallocPt(queryPt_);
    queryPt_ = NULL;
  }
  if (kdTree_fixed_ != NULL) {
    delete kdTree_fixed_;
    kdTree_fixed_ = NULL;
  }
  if (dataPts_fixed_ != NULL) {
    annDeallocPts(dataPts_fixed_);
    dataPts_fixed_ = NULL;
  }
}

void ANNDataStructure::InitializeANN(const PointCloud& pc,
  const double& weightNormal,
  const double& weightColor) {
  dataPts_fixed_ = annAllocPts(static_cast<int> (pc.GetPointArray()->size()), 6);
  for (unsigned v_id = 0; v_id < pc.GetPointArray()->size();
    ++v_id) {
    const Surfel3D& point = (*pc.GetPointArray())[v_id];
    for (int i = 0; i < 3; i++) {
      dataPts_fixed_[v_id][i] = point.position[i];
      dataPts_fixed_[v_id][i + 3] = point.normal[i] * weightNormal;
      dataPts_fixed_[v_id][i + 6] = point.color[i] * weightColor;
    }
  }
  kdTree_fixed_ = new ANNkd_tree(dataPts_fixed_,
    static_cast<int> (pc.GetPointArray()->size()), 6);

  queryPt_ = annAllocPt(9);
  nnIdx_ = new ANNidx[1];
  dists_ = new ANNdist[1];
}

void JointPairwiseRegistration::Compute(const vector<PointCloud>& scans,
  const float& max_distance,
  const float& overlap_ratio,
  const float& weight_normal,
  const float& weight_color,
  const unsigned& down_sampling_rate,
  const unsigned& num_reweighting_iters,
  const unsigned& num_gauss_newton_iters,
  const float& weight_point2planeDis,
  vector<Affine3d>* opt_poses,
  vector<PointCorres2> *pointcorres) {
  Initialize_ANN(scans, weight_normal, weight_color);
  for (unsigned iter = 0; iter < num_reweighting_iters; ++iter) {
    Compute_Correspondences(scans,
      *opt_poses, 
      max_distance,
      overlap_ratio,
      weight_normal,
      weight_color,
      down_sampling_rate,
      pointcorres);
    Correspondence_Reweighting(pointcorres);
    Registration_With_KnownCorres(scans, *pointcorres,
      num_gauss_newton_iters,
      weight_point2planeDis,
      opt_poses);
  }
}
 
void JointPairwiseRegistration::Registration_With_KnownCorres(
  const vector<PointCloud>& scans,
  const vector<PointCorres2>& pointcorres,
  const unsigned& num_gauss_newton_iters,
  const float& weight_point2planeDis,
  vector<Affine3d>* opt_poses) {
  int numscans = static_cast<int> (scans.size());
  int dim = 6 * numscans;
  DMatrixD matA_full(dim, dim);
  DVectorD vecb_full(dim);
  
  // matrix data structure for encoding the linear systems
  DMatrixD matA(dim - 6, dim - 6);
  DVectorD vecb(dim - 6);
  DVectorD vecx(dim - 6);
  double vecJ_plane[12], vecJ_point[3][12], velocity[6];
  unsigned idx_trans[12];
  // Initialize
  for (unsigned i = 0; i < 12; ++i) {
    vecJ_plane[i] = 0.0;
    vecJ_point[0][i] = vecJ_point[1][i] = vecJ_point[2][i] = 0.0;
  }
  for (unsigned i = 0; i < 3; ++i) {
    vecJ_point[i][i] = 1.0;
    vecJ_point[i][6+i] = -1.0;
  }
  // Buffer for encoding the relative transformations
  vector<Affine3d> relativetrans;
  relativetrans.resize(numscans*numscans);

  double w1 = weight_point2planeDis, w2 = 1.0 - weight_point2planeDis;

  for (unsigned iter = 0; iter < num_gauss_newton_iters; ++iter) {
    matA_full.SetZero();
    vecb_full.SetZero();

    //Pre-compute the relative transforms
    for (unsigned sId = 0; sId < numscans; ++sId) {
      const Affine3d& source = (*opt_poses)[sId];
      for (unsigned tId = 0; tId < numscans; ++tId) {
        const Affine3d& target = (*opt_poses)[tId];
        relativetrans[sId*numscans+tId] = target.Inverse()*source;
      }
    }
    double energy = 0.0;
    for (unsigned corId = 0; corId < pointcorres.size(); ++corId) {
      const PointCorres2& cor = pointcorres[corId];
      unsigned offset = cor.sourceSurfId * numscans + cor.targetSurfId;
      for (unsigned i = 0; i < 6; ++i) {
        idx_trans[i] = 6 * cor.sourceSurfId + i;
        idx_trans[i + 6] = 6 * cor.targetSurfId + i;
      }
      const Affine3d& pose = relativetrans[offset];
      const Surfel3D& query =
        (*scans[cor.sourceSurfId].GetPointArray())[cor.sourcePointId];
      const Surfel3D& foot =
        (*scans[cor.targetSurfId].GetPointArray())[cor.targetPointId];
      Vector3d queryPos = pose[0] + pose[1]*query.position[0]
        + pose[2]*query.position[1] + pose[3]*query.position[2];
      Vector3d pointDisVec;
      pointDisVec[0] = queryPos[0] - foot.position[0];
      pointDisVec[1] = queryPos[1] - foot.position[1];
      pointDisVec[2] = queryPos[2] - foot.position[2];
      double planeDis = pointDisVec[0] * query.normal[0]
        + pointDisVec[1] * query.normal[1]
        + pointDisVec[2] * query.normal[2];
      vecJ_plane[0] = foot.normal[0];
      vecJ_plane[1] = foot.normal[1];
      vecJ_plane[2] = foot.normal[2];
      vecJ_plane[3] = queryPos[1] * foot.normal[2] - queryPos[2] * foot.normal[1];
      vecJ_plane[4] = queryPos[2] * foot.normal[0] - queryPos[0] * foot.normal[2];
      vecJ_plane[5] = queryPos[0] * foot.normal[1] - queryPos[1] * foot.normal[0];
      vecJ_point[0][4] = queryPos[2];
      vecJ_point[0][5] = -queryPos[1];
      vecJ_point[1][3] = -queryPos[2];
      vecJ_point[1][5] = queryPos[0];
      vecJ_point[2][3] = queryPos[1];
      vecJ_point[2][4] = -queryPos[0];
      vecJ_point[0][10] = -queryPos[2];
      vecJ_point[0][11] = queryPos[1];
      vecJ_point[1][9] = queryPos[2];
      vecJ_point[1][11] = -queryPos[0];
      vecJ_point[2][9] = -queryPos[1];
      vecJ_point[2][10] = queryPos[0];
      // 
      for (unsigned i = 0; i < 6; ++i)
        vecJ_plane[6 + i] = -vecJ_plane[i];

      // Update the Hessian matrix
      for (unsigned i = 0; i < 12; ++i) {
        for (unsigned j = i; j < 12; ++j) {
          float tp_val = vecJ_point[0][i] * vecJ_point[0][j]
            + vecJ_point[1][i] * vecJ_point[1][j]
            + vecJ_point[2][i] * vecJ_point[2][j];
          matA_full[idx_trans[i]][idx_trans[j]] +=
            (w1 * vecJ_plane[i] * vecJ_plane[j] + w2 * tp_val)*cor.weight;
        }
        vecb_full[idx_trans[i]] -= (w1 * planeDis * vecJ_plane[i]
          + w2 * (pointDisVec[0] * vecJ_point[0][i]
            + pointDisVec[1] * vecJ_point[1][i]
            + pointDisVec[2] * vecJ_point[2][i]))* cor.weight;
      }
      // update the energy
      energy += (w1 * planeDis * planeDis + w2 * pointDisVec.getSqrNorm()) * cor.weight;
    }
    for (int i = 0; i < dim; ++i)
      for (int j = 0; j < i; ++j)
        matA_full[i][j] = matA_full[j][i];
 
    for (int i = 0; i < dim - 6; ++i) {
      vecb[i] = vecb_full[i + 6];
      for (int j = 0; j < dim - 6; ++j) {
        matA[i][j] = matA_full[i + 6][j + 6];
      }
    }
    // Perform line search
    double normOfA = 0.0;
    for (int i = 0; i < dim - 6; ++i)
      normOfA += matA[i][i];
    normOfA /= double(dim - 6.0);
    normOfA = normOfA * 1e-3;

    SolveNxN(matA, vecb, &vecx);
    bool search_success = false;
    for (unsigned searchIter = 0; searchIter < 10; searchIter++) {
      vector<Affine3d> cand_poses = *opt_poses;
      for (unsigned id = 1; id < numscans; ++id) {
        for (unsigned i = 0; i < 6; ++i)
          velocity[i] = vecx[(id - 1) * 6 + i];
        Affine3d pose_diff(velocity);
        cand_poses[id] = pose_diff*cand_poses[id];
        cand_poses[id].Normalize();
      }
      double energy_new = RegistrationObjectiveValue(scans,
        pointcorres,
        cand_poses,
        weight_point2planeDis);
      if (energy_new < energy) {
        *opt_poses = cand_poses;
        search_success = true;
        break;
      } else {
        for (int i = 0; i < dim - 6; ++i)
          matA[i][i] += normOfA;
        normOfA *= 2.0;
        SolveNxN(matA, vecb, &vecx);
      }
    }
  }
}

double JointPairwiseRegistration::RegistrationObjectiveValue(
  const vector<PointCloud>& scans,
  const vector<PointCorres2>& pointcorres,
  const vector<Affine3d>& opt_poses,
  const float& weight_point2planeDis) {
  int numscans = static_cast<unsigned> (scans.size());
  vector<Affine3d> relativetrans;
  relativetrans.resize(numscans*numscans);
  for (unsigned sId = 0; sId < numscans; ++sId) {
    const Affine3d& source = opt_poses[sId];
    for (unsigned tId = 0; tId < numscans; ++tId) {
      const Affine3d& target = opt_poses[tId];
      relativetrans[sId*numscans+tId] = target.Inverse() * source;
    }
  }

  double w1 = weight_point2planeDis,
    w2 = 1.0- weight_point2planeDis,
    energy = 0.0;
  for (unsigned corId = 0; corId < pointcorres.size(); ++corId) {
    const PointCorres2& cor = pointcorres[corId];
    unsigned offset = cor.sourceSurfId * numscans + cor.targetSurfId;
    const Affine3d& pose = relativetrans[offset];
    const Surfel3D& query =
      (*scans[cor.sourceSurfId].GetPointArray())[cor.sourcePointId];
    const Surfel3D& foot =
      (*scans[cor.targetSurfId].GetPointArray())[cor.targetPointId];
    Vector3d queryPos = pose[0] + pose[1] * query.position[0]
      + pose[2] * query.position[1] + pose[3] * query.position[2];
    Vector3d pointDisVec;
    pointDisVec[0] = queryPos[0] - foot.position[0];
    pointDisVec[1] = queryPos[1] - foot.position[1];
    pointDisVec[2] = queryPos[2] - foot.position[2];
    double planeDis = pointDisVec[0] * query.normal[0]
      + pointDisVec[1] * query.normal[1]
      + pointDisVec[2] * query.normal[2];
    // update the energy
    energy += (w1 * planeDis * planeDis 
      + w2 * pointDisVec.getSqrNorm()) * cor.weight;
  }
  return energy;
}

void JointPairwiseRegistration::Initialize_ANN(
  const vector<PointCloud>& scans,
  const float &weightNormal,
  const float &weightColor) {
  for (unsigned id = 0; id < nn_search_dss_.size(); ++id)
    nn_search_dss_[id].ClearANN();
  nn_search_dss_.clear();
  nn_search_dss_.resize(scans.size());
  for (unsigned id = 0; id < nn_search_dss_.size(); ++id)
    nn_search_dss_[id].InitializeANN(scans[id], weightNormal, weightColor);
}

void JointPairwiseRegistration::Clear_ANN() {
  for (unsigned id = 0; id < nn_search_dss_.size(); ++id)
    nn_search_dss_[id].ClearANN();
}


void JointPairwiseRegistration::Correspondence_Reweighting(
  vector<PointCorres2>* pointcorres) {
  vector<float> temp;
  temp.resize(pointcorres->size());
  for (unsigned i = 0; i < pointcorres->size(); ++i)
    temp[i] = (*pointcorres)[i].weight;
  sort(temp.begin(), temp.end());
  float sigma2 = temp[static_cast<unsigned> (temp.size()*0.5)];
  for (unsigned i = 0; i < pointcorres->size(); ++i)
    //    (*pointcorres)[i].weight = exp(-(*pointcorres)[i].weight / 2 / sigma2);
    (*pointcorres)[i].weight = sigma2 / (sigma2 + (*pointcorres)[i].weight);
}

void JointPairwiseRegistration::Compute_Correspondences(
  const vector<PointCloud>& scans,
  const vector<Affine3d> &cur_poses,
  const float& max_distance,
  const float& overlap_ratio,
  const float& weight_normal,
  const float& weight_color,
  const unsigned& down_sampling_rate,
  vector<PointCorres2>* pointcorres) {
  // Generate samples for pair-wise registration based on a down-sampling rate
  vector<vector<unsigned>> sampleIndices;
  sampleIndices.resize(scans.size());
  unsigned total_num_samples = 0;
  for (unsigned scanId = 0; scanId < scans.size(); ++scanId) {
    const PointCloud& pc = scans[scanId];
    vector<unsigned>* indices = &sampleIndices[scanId];
    unsigned num_samples = static_cast<unsigned> (pc.GetPointArray()->size()
      / down_sampling_rate);
    total_num_samples += num_samples;
    indices->resize(num_samples);
    for (unsigned i = 0; i < num_samples; ++i) {
      (*indices)[i] = i * down_sampling_rate;
    }
  }
  unsigned num_scans = static_cast<unsigned> (scans.size());
  unsigned max_num_corres = (num_scans - 1) * total_num_samples;
  unsigned num_corres = 0;
  pointcorres->resize(max_num_corres);
  for (unsigned fixedSurfId = 0; fixedSurfId < scans.size(); ++fixedSurfId) {
    Affine3d fixed_pose_inv = cur_poses[fixedSurfId].Inverse();
    ANNDataStructure* ann_search = &nn_search_dss_[fixedSurfId];
    const PointCloud& fixed_pc = scans[fixedSurfId];
    for (unsigned movingSurfId = 0; movingSurfId < scans.size(); ++movingSurfId) {
      if (fixedSurfId == movingSurfId)
        continue;

      ANNDataStructure* ann_search_moving = &nn_search_dss_[movingSurfId];
      Affine3d cur_pose = fixed_pose_inv * cur_poses[movingSurfId];
      Affine3d cur_pose_inv = cur_pose.Inverse();
      const PointCloud& moving_pc = scans[movingSurfId];
      const vector<unsigned>& indices = sampleIndices[movingSurfId];

      vector<PointCorres2> cand_corres;
      cand_corres.resize(indices.size());
      unsigned num_valid_corres = 0;
      for (unsigned i = 0; i < indices.size(); ++i) {
        const Surfel3D& query = (*moving_pc.GetPointArray())[indices[i]];
        Vector3d pt_pos_transformed = cur_pose[0]
          + cur_pose[1] * query.position[0]
          + cur_pose[2] * query.position[1]
          + cur_pose[3] * query.position[2];
        Vector3d pt_nor_transformed = cur_pose[1] * query.normal[0]
          + cur_pose[2] * query.normal[1]
          + cur_pose[3] * query.normal[2];

        for (int k = 0; k < 3; ++k) {
          ann_search->queryPt_[k] = pt_pos_transformed[k];
          ann_search->queryPt_[k + 3] = weight_normal*pt_nor_transformed[k];
          ann_search->queryPt_[k+6] = weight_color*query.color[k];
        }
        ann_search->kdTree_fixed_->annkSearch(
          ann_search->queryPt_,
          1,
          ann_search->nnIdx_,
          ann_search->dists_,
          0.0);
        const Surfel3D& foot = (*fixed_pc.GetPointArray())[ann_search->nnIdx_[0]];
        //
        for (unsigned k = 0; k < 3; ++k)
            pt_pos_transformed[k] -= foot.position[k];
        double sqrDis = pt_pos_transformed.getSqrNorm();
        //
        if (sqrDis < max_distance * max_distance) {
          cand_corres[num_valid_corres].sourcePointId = indices[i];
          cand_corres[num_valid_corres].sourceSurfId = movingSurfId;
          cand_corres[num_valid_corres].targetPointId = ann_search->nnIdx_[0];
          cand_corres[num_valid_corres].targetSurfId = fixedSurfId;
          cand_corres[num_valid_corres].weight = static_cast<float> (sqrDis);
          num_valid_corres++;
        }
      }
      if (num_valid_corres / double(indices.size()) > overlap_ratio) {
        for (unsigned i = 0; i < num_valid_corres; ++i) {
          (*pointcorres)[num_corres] = cand_corres[i];
          num_corres++;
        }
      }
    }
  }
  pointcorres->resize(num_corres);
}

bool JointPairwiseRegistration::SolveNxN(const DMatrixD& matA, const DVectorD& vecb,
  DVectorD* vecx) {
  // Using LLT factorization to solve the symmetric linear system
  const double fTolerance = 1e-20;
  int dim = vecb.GetDim();
  if (matA.GetColsDim() != dim || matA.GetRowsDim() != dim)
    return false;

  DMatrixD Lower(dim, dim);
  DVectorD afV(dim);
  vecx->SetDim(dim);
  //
  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j)
      Lower[i][j] = matA[i][j];
    (*vecx)[i] = vecb[i];
  }

  for (int i1 = 0; i1 < dim; ++i1) {
    for (int i0 = 0; i0 < i1; ++i0)
      afV[i0] = Lower[i1][i0] * Lower[i0][i0];

    afV[i1] = Lower[i1][i1];
    for (int i0 = 0; i0 < i1; ++i0)
      afV[i1] -= Lower[i1][i0] * afV[i0];

    Lower[i1][i1] = afV[i1];
    if (fabs(afV[i1]) <= fTolerance) //singular
      return false;

    double fInv = 1.0f / afV[i1];
    for (int i0 = i1 + 1; i0 < dim; ++i0) {
      for (int i2 = 0; i2 < i1; ++i2)
        Lower[i0][i1] -= Lower[i0][i2] * afV[i2];
      Lower[i0][i1] *= fInv;
    }
  }

  // Solve Ax = B.
  // Forward substitution
  for (int i0 = 0; i0 < dim; ++i0) {
    for (int i1 = 0; i1 < i0; ++i1)
      (*vecx)[i0] -= Lower[i0][i1] * (*vecx)[i1];
  }

  // Diagonal division:  Let y = L^t x, then Dy = z.  Algorithm stores
  // y terms in B vector.
  for (int i0 = 0; i0 < dim; ++i0) {
    if (fabs(Lower[i0][i0]) <= fTolerance)
      return false;
    (*vecx)[i0] /= Lower[i0][i0];
  }

  // Back substitution:  Solve L^t x = y.  Algorithm stores x terms in
  // B vector.
  for (int i0 = dim-2; i0 >= 0; i0--) {
    for (int i1 = i0 + 1; i1 < dim; ++i1)
      (*vecx)[i0] -= Lower[i1][i0] * (*vecx)[i1];
  }
  return true;
}