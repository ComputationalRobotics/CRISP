// NOTE: the solver uses CSR format for the sparse matrix operation, which is super faster than triplet format, but require more careful handling when use the memory copy.
// The solver rely on a third party sparse convex QP solver to solve the generated subproblem. 
#ifndef SOLVER_INTERFACE_H
#define SOLVER_INTERFACE_H
#include "common/BasicTypes.h"
#include "problem_core/OptimizationProblem.h"
#include "solver_core/SolverParameters.h"
#include <boost/filesystem.hpp>
#include "piqp/piqp.hpp" // qp solver
#include <ctime>

namespace CRISP {
class SolverInterface {
public:
    SolverInterface(OptimizationProblem& problem, SolverParameters& parameters) : problem_(problem), solverParameters_(parameters), initialized_(false) {
        
        }
    void initialize(const vector_t& initial_guess) {
        if (!initialized_) {
            initializeProblem(initial_guess);
            initialized_ = true;
        } else {
            resetProblem(initial_guess);
        }
    }

    void initializeProblem(const vector_t& initial_guess) {
        // Initialize the problem
        problemName_ = problem_.getProblemName();
        variableDim_ = problem_.getVariableDim();
        xIterate_.resize(variableDim_);
        xIterateNext_.resize(variableDim_);
        pTrial_.resize(variableDim_);
        xInitial_ = initial_guess;
        xIterate_ = initial_guess;
        currentIterate_ = 0;
        numEqualityConstraints_ = problem_.getNumEqualityConstraints();
        numInequalityConstraints_ = problem_.getNumInequalityConstraints();
        hasEqualityConstraints_ = numEqualityConstraints_ > 0;
        hasInequalityConstraints_ = numInequalityConstraints_ > 0;
        numConstraints_ = numEqualityConstraints_ + numInequalityConstraints_;
        numObjectives_ = problem_.getNumObjectives();
        offsetV_ = variableDim_;
        offsetW_ = offsetV_ + numEqualityConstraints_;
        offsetT_ = offsetW_ + numEqualityConstraints_;
        totalVars_ = variableDim_ + 2 * numEqualityConstraints_ + numInequalityConstraints_;
        subsolution_.resize(totalVars_); // solution of the subproblem
        numNonZerosObjHess_ = problem_.getNumNonZeroObjHessian();
        numNonZerosEqJac_ = problem_.getNumNonZeroEqJacobian();
        numNonZerosIneqJac_ = problem_.getNumNonZeroIneqJacobian();
        // initialize the subproblem data
        subproblem_ = SubproblemData(totalVars_, numEqualityConstraints_, numInequalityConstraints_, variableDim_);
        subproblem_.H.reserve(numNonZerosObjHess_);
        subproblem_.Aeq.reserve(numNonZerosEqJac_ + 2 * numEqualityConstraints_);
        subproblem_.Aineq.reserve(numNonZerosIneqJac_ + numInequalityConstraints_);
        objHessAugCSR_ = CSRSparseMatrix(totalVars_, numNonZerosObjHess_);
        eqJacAugCSR_ = CSRSparseMatrix(numEqualityConstraints_, numNonZerosEqJac_ + 2 * numEqualityConstraints_);
        ineqJacAugCSR_ = CSRSparseMatrix(numInequalityConstraints_, numNonZerosIneqJac_ + numInequalityConstraints_);
        
        // Initialize the parameters
        maxIterations_ = solverParameters_.getParameters("maxIterations")(0);
        trailTol_ = solverParameters_.getParameters("trailTol")(0);
        trustRegionTol_ = solverParameters_.getParameters("trustRegionTol")(0);
        constraintTol_ = solverParameters_.getParameters("constraintTol")(0);
        trustRegionInitRadius_ = solverParameters_.getParameters("trustRegionInitRadius")(0);
        trustRegionMaxRadius_ = solverParameters_.getParameters("trustRegionMaxRadius")(0);
        mu_ = solverParameters_.getParameters("mu")(0);
        weightedMode_ = solverParameters_.getParameters("WeightedMode")(0);
        weightedTol_ = solverParameters_.getParameters("WeightedTolFactor")(0);
        // initialize mu_matrix_ with diagnal mu
        mu_matrix_.resize(numConstraints_, numConstraints_);
        mu_matrix_.reserve(numConstraints_);
        mu_matrix_.setIdentity();
        mu_matrix_ *= mu_;
        mu_matrix_.makeCompressed();
        ones_eq_.resize(numEqualityConstraints_);
        ones_eq_.setOnes();
        ones_ineq_.resize(numInequalityConstraints_);
        ones_ineq_.setOnes();
        muMax_ = solverParameters_.getParameters("muMax")(0);
        etaLow_ = solverParameters_.getParameters("etaLow")(0);
        etaHigh_ = solverParameters_.getParameters("etaHigh")(0);
        trustRegionRadius_ = trustRegionInitRadius_;
        // prepareStaticTriplet();
        if (solverParameters_.getParameters("verbose")(0) > 0) {
            //     std::cout << std::string(60, '=') << '\n';
            std::cout << std::string(60, '=') << '\n';
            std::cout << center("SOLVER INITIALIZATION", 60) << '\n';
            std::cout << std::string(60, '=') << '\n';
            std::cout << "|                      CRISP                                   |\n";
            std::cout << std::string(60, '-') << '\n';
            std::cout << "| Problem Details:                                              \n";
            std::cout << "|   Name:               " << std::setw(30) << std::left << problemName_ << "|\n";
            std::cout << "|   Variables:          " << std::setw(30) << std::left << variableDim_ << "|\n";
            std::cout << "|   Eq Constraints:     " << std::setw(30) << std::left << numEqualityConstraints_ << "|\n";
            std::cout << "|   Ineq Constraints:   " << std::setw(30) << std::left << numInequalityConstraints_ << "|\n";
            std::cout << "|   Total Constraints:  " << std::setw(30) << std::left << numConstraints_ << "|\n";
            std::cout << std::string(60, '-') << '\n';
            std::cout << "| Solver Parameters:                                            \n";
            std::cout << "|   Max Iterations:     " << std::setw(30) << std::left << maxIterations_ << "|\n";
            std::cout << "|   Trail Tolerance:    " << std::setw(30) << std::left << trailTol_ << "|\n";
            std::cout << "|   Trust Reg Tol:      " << std::setw(30) << std::left << trustRegionTol_ << "|\n";
            std::cout << "|   Const Tol:          " << std::setw(30) << std::left << constraintTol_ << "|\n";
            std::cout << "|   Init Trust Rad:     " << std::setw(30) << std::left << trustRegionInitRadius_ << "|\n";
            std::cout << "|   Max Trust Rad:      " << std::setw(30) << std::left << trustRegionMaxRadius_ << "|\n";
            std::cout << "|   Mu:                 " << std::setw(30) << std::left << mu_ << "|\n";
            std::cout << "|   Mu Max:             " << std::setw(30) << std::left << muMax_ << "|\n";
            std::cout << "|   Eta Low:            " << std::setw(30) << std::left << etaLow_ << "|\n";
            std::cout << "|   Eta High:           " << std::setw(30) << std::left << etaHigh_ << "|\n";
            std::cout << "|   Weighted Mode:      " << std::setw(30) << std::left << weightedMode_ << "|\n";
            std::cout << std::string(60, '=') << '\n';
        }
    }
    // Centering string method
    std::string center(const std::string &text, int width) {
        int padding = (width - text.size()) / 2;
        if (padding > 0) {
            return std::string(padding, ' ') + text + std::string(padding, ' ');
        } else {
            return text; // If the text is too long, just return it
        }
    }
    
    void resetProblem(const vector_t& initial_guess) {
        // problem not change, re-solve the problem with different initial_guess and the solver setting
        // reset the solution
        xInitial_ = initial_guess;
        xIterate_ = initial_guess;
        currentIterate_ = 0;
        // reset the parameters
        maxIterations_ = solverParameters_.getParameters("maxIterations")(0);
        trailTol_ = solverParameters_.getParameters("trailTol")(0);
        trustRegionTol_ = solverParameters_.getParameters("trustRegionTol")(0);
        constraintTol_ = solverParameters_.getParameters("constraintTol")(0);
        trustRegionInitRadius_ = solverParameters_.getParameters("trustRegionInitRadius")(0);
        trustRegionMaxRadius_ = solverParameters_.getParameters("trustRegionMaxRadius")(0);
        mu_ = solverParameters_.getParameters("mu")(0);
        muMax_ = solverParameters_.getParameters("muMax")(0);
        etaLow_ = solverParameters_.getParameters("etaLow")(0);
        etaHigh_ = solverParameters_.getParameters("etaHigh")(0);
        trustRegionRadius_ = trustRegionInitRadius_;
        weightedMode_ = solverParameters_.getParameters("WeightedMode")(0);
        weightedTol_ = solverParameters_.getParameters("WeightedTolFactor")(0);
        secondOrderCorrection_ = solverParameters_.getParameters("secondOrderCorrection")(0);
        mu_matrix_.setIdentity();
        mu_matrix_ *= mu_;
        mu_matrix_.makeCompressed();
        // clear history
        // xHistory_.clear();
        // meritHistory_.clear();
        costHistory_.clear();
        // eqViolationHistory_.clear();
        // ineqViolationHistory_.clear();
        // trustRegionRadiusHistory_.clear();
        
    }

    // set the hyperparameters for the solver, like max iterations, trust region radius, etc
    void setHyperParameters(const std::string& name, const vector_t& params) {
        solverParameters_.setParameters(name, params);
    }

    // for update problem parameters like initial guess, terminal states, for integrating in MPC framework.
    void setProblemParameters(const std::string& name, const vector_t& params) {
        problem_.setParameters(name, params);
    }

    // void prepareStaticTriplet() {
    //     // prepare the triplet for the static part of the subproblem
    //     // equality constraints
    //     for (size_t i = 0; i < numEqualityConstraints_; ++i) {
    //         tripletListFixed_Aeq_.emplace_back(i, i + offsetV_, -1.0);
    //         tripletListFixed_Aeq_.emplace_back(i, i + offsetW_, 1.0);
    //     }
    //     // inequality constraints
    //     for (size_t i = 0; i < numInequalityConstraints_; ++i) {
    //         tripletListFixed_Aineq_.emplace_back(i, i + offsetT_, 1.0);
    //     }
    // }



    // print the result

    void solve() {
        // initialization
        obj_ = problem_.evaluateObjective(xIterate_);
        objJac_ = problem_.evaluateObjectiveGradient(xIterate_).toDense().row(0).transpose();
        objHessCSR_ = problem_.evaluateObjectiveHessianCSR(xIterate_);
        objHessMat_.resize(variableDim_, variableDim_);
        objHessMat_.reserve(numNonZerosObjHess_);
        objHessCSR_.toEigenSparseMatrix(objHessMat_);

        // value of the constraints
        eqValues_ = problem_.evaluateEqualityConstraints(xIterate_);
        ineqValues_ = problem_.evaluateInequalityConstraints(xIterate_);
        eqJacCSR_ = problem_.evaluateEqualityConstraintsJacobianCSR(xIterate_);
        eqJacMat_.resize(numEqualityConstraints_, variableDim_);
        eqJacMat_.reserve(numNonZerosEqJac_);
        eqJacCSR_.toEigenSparseMatrix(eqJacMat_);
        ineqJacCSR_ = problem_.evaluateInequalityConstraintsJacobianCSR(xIterate_);
        ineqJacMat_.resize(numInequalityConstraints_, variableDim_);
        ineqJacMat_.reserve(numNonZerosIneqJac_);
        ineqJacCSR_.toEigenSparseMatrix(ineqJacMat_);
        secondOrderCorrectionCount = 0;
        phi_ = evaluateMeritFunction(obj_, eqValues_, ineqValues_);
        q_mu_0_ = evaluateQuadraticModel(obj_, objJac_, eqValues_, ineqValues_, objHessMat_, eqJacMat_, ineqJacMat_);
        time_qp = 0.0;
        time_total = 0.0;
        // main loop, reuse data from the previous iteration to improve efficiency.
        auto startsolve = std::chrono::high_resolution_clock::now();
        for (currentIterate_ = 0; currentIterate_ < maxIterations_; ++currentIterate_) {
            // store the previous iterate
            costHistory_.push_back(obj_);
            if (solverParameters_.getParameters("verbose")(0) > 0) {
                std::cout << "Iteration: " << currentIterate_ << " Objective: " << obj_ << " Merit: " << phi_ << " Trust region: " << trustRegionRadius_ << std::endl;
                std::cout << "Equality violation: " << eqValues_.array().abs().maxCoeff() << " Inequality violation: " << (-ineqValues_).array().maxCoeff() << std::endl;
            }
            // construct the subproblem
            constructSubproblem(objJac_, objHessCSR_, eqValues_, ineqValues_, eqJacCSR_, ineqJacCSR_);
            // auto end = std::chrono::high_resolution_clock::now();
            // std::cout << "Subproblem construction time: " << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() << "us" << std::endl;
            subsolution_ = solveSubproblem(subproblem_); // solve the subproblem
            std::memcpy(pTrial_.data(), subsolution_.data(), variableDim_ * sizeof(scalar_t));
            // evaluate necessary value at the trial step
            xIterateNext_ = xIterate_ + pTrial_;
            // auto starteval = std::chrono::high_resolution_clock::now();
            auto objNext = problem_.evaluateObjective(xIterateNext_);
            auto eqValuesNext = problem_.evaluateEqualityConstraints(xIterateNext_);
            auto ineqValuesNext = problem_.evaluateInequalityConstraints(xIterateNext_);
            // auto endeval = std::chrono::high_resolution_clock::now();
            // std::cout << "Evaluation time: " << std::chrono::duration_cast<std::chrono::microseconds>(endeval - starteval).count() << "us" << std::endl;
            phi_pk_ = evaluateMeritFunction(objNext, eqValuesNext, ineqValuesNext); // mertit function at the trial step
            q_mu_pk_ = evaluateQuadraticModel(obj_, objJac_, eqValues_, ineqValues_, objHessMat_, eqJacMat_, ineqJacMat_, pTrial_); // quadratic model at the trial step
            // second order correction if actual reduction less than 0;
            if (phi_ - phi_pk_ < 0) {
                // std::cout << "actual reduction before second order correction: " << phi_ - phi_pk_ << std::endl;
                // modify the subproblem, resolve for a new trial step to consider the second order correction
                secondOrderCorrectionCount++;
                // change subproblem_.beq->-(eqValuesNext-eqconstraintJac*ptrial) and subproblem_.bineq->-(IneqValuesNext-IneqconstraintJac*ptrial) and resolve the problem.
                subproblem_.beq = -(eqValuesNext - eqJacMat_ * pTrial_);
                subproblem_.bineq = -(ineqValuesNext - ineqJacMat_ * pTrial_);
                subsolution_ = solveSubproblem(subproblem_);
                std::memcpy(pTrial_.data(), subsolution_.data(), variableDim_ * sizeof(scalar_t));
                xIterateNext_ = xIterate_ + pTrial_;
                objNext = problem_.evaluateObjective(xIterateNext_);
                eqValuesNext = problem_.evaluateEqualityConstraints(xIterateNext_);
                ineqValuesNext = problem_.evaluateInequalityConstraints(xIterateNext_);
                q_mu_0_ = evaluateQuadraticModel(obj_, objJac_, -subproblem_.beq, -subproblem_.bineq, objHessMat_, eqJacMat_, ineqJacMat_);
                q_mu_pk_ = evaluateQuadraticModel(obj_, objJac_, -subproblem_.beq, -subproblem_.bineq, objHessMat_, eqJacMat_, ineqJacMat_, pTrial_);
                phi_pk_ = evaluateMeritFunction(objNext, eqValuesNext, ineqValuesNext); // mertit function at the trial step
            }
 
            reduction_ratio_ = (phi_ - phi_pk_) / (q_mu_0_ - q_mu_pk_);
            scalar_t reduction_actual_ = phi_ - phi_pk_;
            if (updateTrustRegionRadius(reduction_ratio_, reduction_actual_)) {
                // if the trial step is accepted, update the iterate
                xIterate_ = xIterateNext_;
                obj_ = objNext;
                eqValues_ = eqValuesNext;
                ineqValues_ = ineqValuesNext;
                phi_ = phi_pk_;
                q_mu_0_ = phi_;
                objJac_ = problem_.evaluateObjectiveGradient(xIterate_).toDense().row(0).transpose();
                objHessCSR_ = problem_.evaluateObjectiveHessianCSR(xIterate_);
                objHessCSR_.toEigenSparseMatrix(objHessMat_);
                eqJacCSR_ = problem_.evaluateEqualityConstraintsJacobianCSR(xIterate_);
                eqJacCSR_.toEigenSparseMatrix(eqJacMat_);
                ineqJacCSR_ = problem_.evaluateInequalityConstraintsJacobianCSR(xIterate_);
                ineqJacCSR_.toEigenSparseMatrix(ineqJacMat_);
            }
            else {
                // if the trial step is rejected, the trust region radius is shrinked but not update the iterate.
                // std::cout<<"Trial step rejected."<<std::endl;
            }

            // check the stopping criteria, increase the penalty if necessary.
            if (checkStoppingCriteria()) {
                break;
            }
        }
        auto endsolve = std::chrono::high_resolution_clock::now();
        time_total = std::chrono::duration_cast<std::chrono::milliseconds>(endsolve - startsolve).count();
        // xHistory_.push_back(xIterate_); // Maintain full state for x
        // meritHistory_.push_back(phi_);
        costHistory_.push_back(obj_);
        // eqViolationHistory_.push_back(eqValues_.array().abs().maxCoeff());
        // ineqViolationHistory_.push_back((-ineqValues_).array().maxCoeff());
        // trustRegionRadiusHistory_.push_back(trustRegionRadius_);
        // saveResults(); // save the results to .mat file
    }

    vector_t getSolution() {
        std::cout << "Optimization problem:" << problemName_ << " solved in " << currentIterate_ << " iterations." << std::endl;
        std::cout << "Solver time: " << time_total << "ms" << std::endl;
        std::cout << "QP solver time: " << time_qp/1000 << "ms" << std::endl;
        std::cout << "max constraint violation: equality: " << eqValues_.array().abs().maxCoeff() << " inequality: " << (-ineqValues_).array().maxCoeff() << std::endl;
        std::cout << "obj from " << costHistory_[0] << " to " << costHistory_.back() << std::endl;
        return xIterate_;
    }

    void saveResults(const std::string& folderPrefix) {
        // // Get the current time as a unique identifier
        // std::time_t t = std::time(nullptr);
        // std::ostringstream fileNameStream;
        // fileNameStream << problemName_ << "_result_" << t;  // Use current time for unique naming
        // std::string folderName = fileNameStream.str();
        // std::string folderPath = folderPrefix + "/" + folderName;

        // try {
        //     // Create necessary directories using Boost
        //     if (!boost::filesystem::exists(folderPrefix)) {
        //         boost::filesystem::create_directory(folderPrefix);
        //     }

        //     boost::filesystem::create_directory(folderPath);

        //     std::string fullFileName = folderPath + "/data.mat";
        //     // MAT-file creation with MATLAB API
        //     MATFile *matFile = matOpen(fullFileName.c_str(), "w");
        //     if (!matFile) {
        //         std::cerr << "Error opening MAT-file for writing: " << fullFileName << std::endl;
        //         return;
        //     }

        //     if (!xHistory_.empty()) {
        //         size_t numVars = xHistory_[0].size();
        //         mxArray* xHistArray = mxCreateDoubleMatrix(numVars, xHistory_.size(), mxREAL);
        //         auto pr = mxGetPr(xHistArray);
        //         for (size_t i = 0; i < xHistory_.size(); ++i) {
        //             std::copy(xHistory_[i].data(), xHistory_[i].data() + numVars, pr + numVars * i);
        //         }
        //         matPutVariable(matFile, "xHistory", xHistArray);
        //         mxDestroyArray(xHistArray);
        //     }
        //     // Helper lambda for writing scalar vectors to MAT-file
        //     auto writeScalarVector = [&](const std::string& varName, const std::vector<scalar_t>& data) {
        //         mxArray* array = mxCreateDoubleMatrix(data.size(), 1, mxREAL);
        //         std::copy(data.begin(), data.end(), mxGetPr(array));
        //         matPutVariable(matFile, varName.c_str(), array);
        //         mxDestroyArray(array);
        //     };

        //     // Save the necessary data
        //     // writeScalarVector("meritHistory", meritHistory_);
        //     writeScalarVector("costHistory", costHistory_);
        //     // writeScalarVector("eqViolationHistory", eqViolationHistory_);
        //     // writeScalarVector("ineqViolationHistory", ineqViolationHistory_);
        //     // writeScalarVector("trustRegionRadiusHistory", trustRegionRadiusHistory_);

        //     matClose(matFile);
        //     } catch (const boost::filesystem::filesystem_error& e) {
        //         std::cerr << "Filesystem error: " << e.what() << std::endl;
        //     }
        //     std::cout << "Results saved to: " << folderPath << std::endl;
    }

private:
    // standard subproblem format for the QP solver.
    struct SubproblemData {
        vector_t g;
        sparse_matrix_t H;
        sparse_matrix_t Aeq;
        sparse_matrix_t Aineq;
        vector_t beq;
        vector_t bineq;
        vector_t lb;
        vector_t ub;
        vector_t x0;
        SubproblemData() = default;
        SubproblemData(size_t totalVars, size_t numEqualityConstraints, size_t numInequalityConstraints, size_t variableDim)
            : g(totalVars),
              H(totalVars, totalVars),
              Aeq(numEqualityConstraints, totalVars),
              Aineq(numInequalityConstraints, totalVars),
              beq(numEqualityConstraints),
              bineq(numInequalityConstraints),
              lb(totalVars),
              ub(totalVars),
              x0(totalVars) 
              {}
    };


    // for piqp
    void constructSubproblem(const vector_t& objJac, const CSRSparseMatrix& objHess, const vector_t& eqValues, const vector_t& ineqValues, const CSRSparseMatrix& eqJac, const CSRSparseMatrix& ineqJac) {
        // construct the subproblem for the QP solver
        // build objecitve gradient and hessian
        if (weightedMode_ > 0){
            subproblem_.g.setOnes();
            subproblem_.g.head(variableDim_) = objJac;
            // use mu_matrix_ to construct g:[objJac;mu_matrix(1:numEqualityConstraints,1:numEqualityConstraints)*ones(numEqualityConstraints,1);mu_matrix(1:numEqualityConstraints,1:numEqualityConstraints)*ones(numEqualityConstraints,1);mu_matrix(numEq:end,numEq:end)*ones(numInequalityConstraints,1)]
            vector_t mu_eq(numEqualityConstraints_,1);
            vector_t mu_ineq(numInequalityConstraints_,1);
            for (size_t i = 0; i < numEqualityConstraints_; ++i) {
                mu_eq[i] = mu_matrix_.valuePtr()[i];
            }
            for (size_t i = 0; i < numInequalityConstraints_; ++i) {
                mu_ineq[i] = mu_matrix_.valuePtr()[numEqualityConstraints_ + i];
            }
            // memory copy
            std::memcpy(subproblem_.g.data() + variableDim_, mu_eq.data(), numEqualityConstraints_ * sizeof(scalar_t));
            std::memcpy(subproblem_.g.data() + offsetW_, mu_eq.data(), numEqualityConstraints_ * sizeof(scalar_t));
            std::memcpy(subproblem_.g.data() + offsetT_, mu_ineq.data(), numInequalityConstraints_ * sizeof(scalar_t));
        }
        else {
            subproblem_.g.setConstant(mu_);
            subproblem_.g.head(variableDim_) = objJac;
        }

        buildBlockObjHessCSR(objHess, objHessAugCSR_);
        objHessAugCSR_.toEigenSparseMatrix(subproblem_.H);
        // build equality constraints
        buildBlockAeqCSR(eqJac, eqJacAugCSR_);
        eqJacAugCSR_.toEigenSparseMatrix(subproblem_.Aeq);
        subproblem_.beq = -eqValues;
        // build inequality constraints
        buildBlockAinCSR(ineqJac, ineqJacAugCSR_);
        ineqJacAugCSR_.toEigenSparseMatrix(subproblem_.Aineq);
        subproblem_.bineq = -ineqValues;
        // bounds
        subproblem_.lb.setZero();
        subproblem_.lb.head(variableDim_).setConstant(-trustRegionRadius_);
        // inf ub bounds for the slack variables
        subproblem_.ub.setConstant(std::numeric_limits<scalar_t>::infinity());
        subproblem_.ub.head(variableDim_).setConstant(trustRegionRadius_);
        subproblem_.x0.setZero();
    }
 

    scalar_t evaluateMeritFunction(const scalar_t& obj, const vector_t& eqValues, const vector_t& ineqValues)
    {
        // return obj + mu_ * (eqValues.array().abs().sum() + ((-ineqValues).array().max(0.0)).sum());
        // use mu_matrix_ to weight the constraints violation
        if (weightedMode_ > 0){
            return obj + ones_eq_.transpose() * mu_matrix_.block(0,0,numEqualityConstraints_,numEqualityConstraints_) * eqValues.array().abs().matrix() + ones_ineq_.transpose() * mu_matrix_.block(numEqualityConstraints_,numEqualityConstraints_,numInequalityConstraints_,numInequalityConstraints_) * (-ineqValues).array().max(0.0).matrix();
        }
        else {
            return obj + mu_ * (eqValues.array().abs().sum() + ((-ineqValues).array().max(0.0)).sum());
        }
    }

    scalar_t evaluateQuadraticModel(const scalar_t& obj, const vector_t& objJac, const vector_t& eqValues, const vector_t& ineqValues, const sparse_matrix_t& objHessMat, const sparse_matrix_t& eqJacMat, const sparse_matrix_t& ineqJacMat, const vector_t& p)
    {
        // return obj + objJac.dot(p) + 0.5 * p.transpose() * objHessMat * p  +  mu_ * ((eqValues + eqJacMat * p).array().abs().sum() + ((-ineqValues - ineqJacMat * p).array().max(0.0)).sum());
        // use mu_matrix_ to weight the constraints violation
        if (weightedMode_ > 0){
            return obj + objJac.dot(p) + 0.5 * p.transpose() * objHessMat * p + ones_eq_.transpose() * mu_matrix_.block(0,0,numEqualityConstraints_,numEqualityConstraints_) * (eqValues + eqJacMat * p).array().abs().matrix() + ones_ineq_.transpose() * mu_matrix_.block(numEqualityConstraints_,numEqualityConstraints_,numInequalityConstraints_,numInequalityConstraints_) * (-ineqValues - ineqJacMat * p).array().max(0.0).matrix();
        }
        else {
            return obj + objJac.dot(p) + 0.5 * p.transpose() * objHessMat * p  +  mu_ * ((eqValues + eqJacMat * p).array().abs().sum() + ((-ineqValues - ineqJacMat * p).array().max(0.0)).sum());
        }
    }

    scalar_t evaluateQuadraticModel(const scalar_t& obj, const vector_t& objJac, const vector_t& eqValues, const vector_t& ineqValues, const sparse_matrix_t& objHessMat, const sparse_matrix_t& eqJacMat, const sparse_matrix_t& ineqJacMat)
    {   
        if (weightedMode_ > 0){
            return obj + ones_eq_.transpose() * mu_matrix_.block(0,0,numEqualityConstraints_,numEqualityConstraints_) * eqValues.array().abs().matrix() + ones_ineq_.transpose() * mu_matrix_.block(numEqualityConstraints_,numEqualityConstraints_,numInequalityConstraints_,numInequalityConstraints_) * (-ineqValues).array().max(0.0).matrix();
            }
        else {
            return obj + mu_ * (eqValues.array().abs().sum() + ((-ineqValues).array().max(0.0)).sum());
        
        }
    }

    // convex QP solver: now using the piqp.
    vector_t solveSubproblem(const SubproblemData& subproblem) {
        auto startsol = std::chrono::high_resolution_clock::now();
        // solve the subproblem using the QP solver
        if (currentIterate_ == 0) {
            // piqpSolver_.settings().max_iter = 200;
            piqpSolver_.setup(subproblem.H, subproblem.g, subproblem.Aeq, subproblem.beq, subproblem.Aineq, subproblem.bineq,
              piqp::nullopt, subproblem.lb, subproblem.ub);
          } else {
            piqpSolver_.update(subproblem.H, subproblem.g, subproblem.Aeq, subproblem.beq, subproblem.Aineq, subproblem.bineq,
              piqp::nullopt, subproblem.lb, subproblem.ub);
          }
        piqp::Status status = piqpSolver_.solve();     
        auto endsol = std::chrono::high_resolution_clock::now();
        time_qp += std::chrono::duration_cast<std::chrono::microseconds>(endsol - startsol).count();
        return piqpSolver_.result().x;

    }

    bool updateTrustRegionRadius(const scalar_t& reduction_ratio, const scalar_t& reduction_actual) {
        if (reduction_ratio < 0 || reduction_actual < 1e-9) {
            trustRegionRadius_ *= 0.25; // Shrink trust region
            // std::cout << "Step not accepted, shrinking trust region.\n";
            return false; // Step not accepted
        } else if (reduction_ratio < etaLow_) {
            trustRegionRadius_ *= 0.25; // shrink
        } else if (reduction_ratio > etaHigh_ && pTrial_.array().abs().maxCoeff() > 0.95 * trustRegionRadius_) {
            // std::cout << "Increasing the trust region.\n";
            trustRegionRadius_ = std::min(2 * trustRegionRadius_, trustRegionMaxRadius_); // Increase trust region
        } else {
            // std::cout << "Trust region unchanged.\n";
        }
        return true; // Step accepted

    }

    bool checkStoppingCriteria() {
        // check the stopping criteria
        if (trustRegionRadius_ < trustRegionTol_ || pTrial_.norm()/xIterate_.norm() < trailTol_) {
            std::cout << "current merit function converge, examing the constraints violation.\n" << std::endl;
            if (eqValues_.array().abs().maxCoeff() < constraintTol_ && (-ineqValues_).array().maxCoeff() < constraintTol_) {
                std::cout << "Optimization converged.\n" << std::endl;
                return true;
            }
            else {
                if (weightedMode_ > 0){
                    scalar_t max_mu = std::numeric_limits<scalar_t>::lowest();

                    // Iterate over diagonal elements only
                    for (int i = 0; i < numConstraints_; ++i) {
                        if (mu_matrix_.valuePtr()[i] > max_mu) {
                            max_mu = mu_matrix_.valuePtr()[i];
                        }
                    }

                    if (max_mu == muMax_) {
                        std::cout << "penalty maxed out, check the solution.\n" << std::endl;
                        return true;
                    }
                    // mu_ = std::min(10 * mu_, muMax_);
                    // update mu_matrix_ according to which constraints are violated.
                    for (size_t i = 0; i < numEqualityConstraints_; ++i) {
                        if (eqValues_.array().abs()[i] > weightedTol_*constraintTol_) {
                            mu_matrix_.valuePtr()[i] = std::min(10 * mu_matrix_.valuePtr()[i], muMax_);
                            // std::cout << "equality constraint " << i << " violated, increase penalty to " << mu_matrix_.valuePtr()[i] << std::endl;
                
                        }
                    }
                    for (size_t i = 0; i < numInequalityConstraints_; ++i) {
                        if (-ineqValues_[i] > weightedTol_*constraintTol_) {
                            mu_matrix_.valuePtr()[numEqualityConstraints_ + i] = std::min(10 * mu_matrix_.valuePtr()[numEqualityConstraints_ + i], muMax_);
                            // std::cout << "inequality constraint " << i << " violated, increase penalty to " << mu_matrix_.valuePtr()[numEqualityConstraints_ + i] << std::endl;
                        }
                    }
                }
                else {
                    if (mu_ == muMax_) {
                        std::cout << "penalty maxed out, check the solution.\n" << std::endl;
                        return true;
                    }
                    // std::cout << "increase penalty" << std::endl;
                    mu_ = std::min(10 * mu_, muMax_);
                }
                phi_ = evaluateMeritFunction(obj_, eqValues_, ineqValues_);
                q_mu_0_ = evaluateQuadraticModel(obj_, objJac_, eqValues_, ineqValues_, objHessMat_, eqJacMat_, ineqJacMat_);
                // trustRegionRadius_ = trustRegionInitRadius_;
            }
        }
        return false;
    }

    // save the results to a matlab file


    // [A:offsetV:-I:offsetW:I:offsetT:0]: needs careful handling.
    void buildBlockAeqCSR(const CSRSparseMatrix& Aeq, CSRSparseMatrix& Aeq_aug){
        size_t numNonZeroCurrentRow;
        size_t numNonZeroTotal = 0;
        for (size_t i = 0; i < numEqualityConstraints_; ++i) {
            Aeq_aug.outerIndex[i+1] = Aeq.outerIndex[i+1] + 2 * (i+1); // outer iterator
            numNonZeroCurrentRow = Aeq.outerIndex[i + 1] - Aeq.outerIndex[i];
            // Aeq_aug[1].segment(numNonZeroTotal, numNonZeroCurrentRow + 2) << Aeq[1].segment(Aeq[0](i), numNonZeroCurrentRow), i + offsetV_, i + offsetW_;
            std::memcpy(Aeq_aug.innerIndices.data() + numNonZeroTotal, Aeq.innerIndices.data() + Aeq.outerIndex[i], numNonZeroCurrentRow * sizeof(size_t));
            // add i + offsetV_ and i + offsetW_ to the end of the vector
            Aeq_aug.innerIndices[numNonZeroTotal + numNonZeroCurrentRow] = i + offsetV_;
            Aeq_aug.innerIndices[numNonZeroTotal + numNonZeroCurrentRow + 1] = i + offsetW_;
            // Aeq_aug[2].segment(numNonZeroTotal, numNonZeroCurrentRow + 2) << Aeq[2].segment(Aeq[0](i), numNonZeroCurrentRow), -1.0, 1.0;
            std::memcpy(Aeq_aug.values.data() + numNonZeroTotal, Aeq.values.data() + Aeq.outerIndex[i], numNonZeroCurrentRow * sizeof(scalar_t));
            Aeq_aug.values[numNonZeroTotal + numNonZeroCurrentRow] = -1.0;
            Aeq_aug.values[numNonZeroTotal + numNonZeroCurrentRow + 1] = 1.0;
            numNonZeroTotal += numNonZeroCurrentRow + 2;
        }
    }
    // [A:offsetT:I]: needs careful handling.
    void buildBlockAinCSR(const CSRSparseMatrix& Aineq, CSRSparseMatrix& Aineq_aug){
        size_t numNonZeroCurrentRow;
        size_t numNonZeroTotal = 0;
        for (size_t i = 0; i < numInequalityConstraints_; ++i) {
            Aineq_aug.outerIndex[i+1] = Aineq.outerIndex[i+1] + (i+1); // outer iterator
            numNonZeroCurrentRow = Aineq.outerIndex[i + 1] - Aineq.outerIndex[i];
            // Aineq_aug[1].segment(numNonZeroTotal, numNonZeroCurrentRow + 1) << Aineq[1].segment(Aineq[0](i), numNonZeroCurrentRow), i + offsetT_;
            std::memcpy(Aineq_aug.innerIndices.data() + numNonZeroTotal, Aineq.innerIndices.data() + Aineq.outerIndex[i], numNonZeroCurrentRow * sizeof(size_t));
            Aineq_aug.innerIndices[numNonZeroTotal + numNonZeroCurrentRow] = i + offsetT_;
            // Aineq_aug[2].segment(numNonZeroTotal, numNonZeroCurrentRow + 1) << Aineq[2].segment(Aineq[0](i), numNonZeroCurrentRow), 1.0;
            std::memcpy(Aineq_aug.values.data() + numNonZeroTotal, Aineq.values.data() + Aineq.outerIndex[i], numNonZeroCurrentRow * sizeof(scalar_t));
            Aineq_aug.values[numNonZeroTotal + numNonZeroCurrentRow] = 1.0;
            numNonZeroTotal += numNonZeroCurrentRow + 1;
        }
    }
    // [H,0;0,0]
    void buildBlockObjHessCSR(const CSRSparseMatrix& objHess, CSRSparseMatrix& objHess_aug){
        std::memcpy(objHess_aug.innerIndices.data(), objHess.innerIndices.data(), objHess.innerIndices.size() * sizeof(size_t));
        std::memcpy(objHess_aug.values.data(), objHess.values.data(), objHess.values.size() * sizeof(scalar_t));
        // objHess_aug[0].segment(0, variableDim_ + 1) = objHess[0].segment(0, variableDim_ + 1);
        // objHess_aug[0].segment(variableDim_ + 1, 2 * numEqualityConstraints_ + numInequalityConstraints_) = vector_t::Constant(2 * numEqualityConstraints_ + numInequalityConstraints_, objHess[0](variableDim_));
        std::memcpy(objHess_aug.outerIndex.data(), objHess.outerIndex.data(), (variableDim_ + 1) * sizeof(size_t));
        SizeVector temp(2 * numEqualityConstraints_ + numInequalityConstraints_, objHess.outerIndex[variableDim_]);
        std::memcpy(objHess_aug.outerIndex.data() + variableDim_ + 1, temp.data(), temp.size() * sizeof(size_t));
    }

    // ----- variables ----- //
    std::string problemName_;
    piqp::SparseSolver<scalar_t> piqpSolver_;
    SubproblemData subproblem_;
    OptimizationProblem problem_;
    SolverParameters solverParameters_;
    bool initialized_;
    size_t variableDim_;
    size_t secondOrderCorrectionCount;
    vector_t subsolution_;
    vector_t xIterate_;
    vector_t xIterateTemp_;
    vector_t xInitial_;
    vector_t pTrial_;
    vector_t xIterateNext_;
    vector_t eqValues_;
    vector_t ineqValues_;
    vector_t objJac_;
    vector_t ones_eq_;
    vector_t ones_ineq_;
    // triplet_vector_t objHessTriplets_;
    // triplet_vector_t eqJacTriplets_;
    // triplet_vector_t ineqJacTriplets_;
    sparse_matrix_t objHessMat_;
    sparse_matrix_t eqJacMat_;
    sparse_matrix_t ineqJacMat_;
    sparse_matrix_t mu_matrix_;
    CSRSparseMatrix objHessCSR_;
    CSRSparseMatrix eqJacCSR_;
    CSRSparseMatrix ineqJacCSR_;
    CSRSparseMatrix objHessAugCSR_;
    CSRSparseMatrix eqJacAugCSR_;
    CSRSparseMatrix ineqJacAugCSR_;
    size_t currentIterate_;
    size_t numEqualityConstraints_;
    size_t numInequalityConstraints_;
    size_t numConstraints_;
    size_t numObjectives_;
    size_t maxIterations_;
    size_t numNonZerosObjHess_;
    size_t numNonZerosEqJac_;
    size_t numNonZerosIneqJac_;
    scalar_t obj_; //objective function
    scalar_t phi_; //merit function
    scalar_t phi_pk_; //merit function at next step
    scalar_t q_mu_0_; //quadratic model 
    scalar_t q_mu_pk_; //estimation at next step
    scalar_t reduction_ratio_; //reduction ratio between merit function and the quadratic model.
    scalar_t trailTol_;
    scalar_t trustRegionTol_;
    scalar_t constraintTol_;
    scalar_t trustRegionInitRadius_;
    scalar_t trustRegionMaxRadius_;
    scalar_t mu_;
    scalar_t muMax_;
    scalar_t etaLow_;
    scalar_t etaHigh_;
    scalar_t trustRegionRadius_;
    scalar_t time_qp;
    scalar_t time_total;
    scalar_t initial_squareNorm_eq_;
    scalar_t initial_squareNorm_ineq_;
    scalar_t weightedMode_;
    scalar_t secondOrderCorrection_;
    scalar_t weightedTol_;
    size_t offsetV_; //subproblem: offset for the slack variables for the equality constraints.
    size_t offsetW_; //subproblem: offset for the slack variables for the inequality constraints.
    size_t offsetT_; //subproblem: offset for the slack variables for the inequality constraints.
    size_t totalVars_; //total number of variables in the subproblem.
    size_t increaseCount_; 
    bool hasEqualityConstraints_;
    bool hasInequalityConstraints_;
    triplet_vector_t tripletListFixed_Aeq_;
    triplet_vector_t tripletListFixed_Aineq_;
    // // history
    // std::vector<vector_t> xHistory_;
    // std::vector<scalar_t> meritHistory_;
    std::vector<scalar_t> costHistory_;
    // std::vector<scalar_t> eqViolationHistory_;
    // std::vector<scalar_t> ineqViolationHistory_;
    // std::vector<scalar_t> trustRegionRadiusHistory_;
};
} // namespace CRISP
#endif // SOLVER_INTERFACE_H