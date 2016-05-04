## Process this file with automake to produce Makefile.in

## Using wildcards in these lists doesn't work. Also GNU make's
## ($wildcard,) doesn't work. It gives warning message about
## portability, but in the end doesn't work, I tried :-).
REGFILES = main.cpp					\
 invsigma.h invsigma.cpp				\
 gendata.h gendata.cpp					\
 eigen_mematrix.h eigen_mematrix.cpp			\
 command_line_settings.h command_line_settings.cpp	\
 reg1.h reg1.cpp					\
 usage.h usage.cpp					\
 utilities.h utilities.cpp				\
 phedata.h phedata.cpp					\
 regdata.h regdata.cpp					\
 maskedmatrix.h maskedmatrix.cpp			\
 main_functions_dump.h main_functions_dump.cpp		\
 mlinfo.h mlinfo.cpp


## Filevector files:
FVSRC = fvlib/AbstractMatrix.cpp fvlib/CastUtils.cpp			\
 fvlib/convert_util.cpp fvlib/FileVector.cpp fvlib/FilteredMatrix.cpp	\
 fvlib/frutil.cpp fvlib/Logger.cpp fvlib/RealHandlerWrapper.cpp		\
 fvlib/ReusableFileHandle.cpp fvlib/Transposer.cpp
FVHEADERS = fvlib/AbstractMatrix.h fvlib/CastUtils.h fvlib/const.h	\
 fvlib/convert_util.h fvlib/FileVector.h fvlib/frutil.h			\
 fvlib/frversion.h fvlib/Logger.h fvlib/Transposer.h			\
 fvlib/ReusableFileHandle.h fvlib/RealHandlerWrapper.h			\
 fvlib/FilteredMatrix.h

## R headers, used by pacoxph
RHEADERS = include/R.h include/Rmath.h include/R_ext/Arith.h		\
 include/R_ext/Boolean.h include/Rconfig.h include/R_ext/Complex.h	\
 include/R_ext/Constants.h include/R_ext/Error.h			\
 include/R_ext/libextern.h include/R_ext/Memory.h			\
 include/R_ext/Print.h include/R_ext/Random.h include/R_ext/Utils.h	\
 include/R_ext/RS.h

## Eigen files, listed here so they get included in the source tar-ball.
EIGENFILES =									\
 eigen-3.2.1/Eigen/Jacobi							\
 eigen-3.2.1/Eigen/IterativeLinearSolvers					\
 eigen-3.2.1/Eigen/MetisSupport							\
 eigen-3.2.1/Eigen/Array							\
 eigen-3.2.1/Eigen/SparseLU							\
 eigen-3.2.1/Eigen/LU								\
 eigen-3.2.1/Eigen/SparseCholesky						\
 eigen-3.2.1/Eigen/Core								\
 eigen-3.2.1/Eigen/SPQRSupport							\
 eigen-3.2.1/Eigen/UmfPackSupport						\
 eigen-3.2.1/Eigen/Geometry							\
 eigen-3.2.1/Eigen/StdVector							\
 eigen-3.2.1/Eigen/Sparse							\
 eigen-3.2.1/Eigen/OrderingMethods						\
 eigen-3.2.1/Eigen/StdDeque							\
 eigen-3.2.1/Eigen/Eigen2Support						\
 eigen-3.2.1/Eigen/Householder							\
 eigen-3.2.1/Eigen/Eigenvalues							\
 eigen-3.2.1/Eigen/SparseCore							\
 eigen-3.2.1/Eigen/QtAlignedMalloc						\
 eigen-3.2.1/Eigen/PaStiXSupport						\
 eigen-3.2.1/Eigen/StdList							\
 eigen-3.2.1/Eigen/Dense							\
 eigen-3.2.1/Eigen/PardisoSupport						\
 eigen-3.2.1/Eigen/SparseQR							\
 eigen-3.2.1/Eigen/SVD								\
 eigen-3.2.1/Eigen/CMakeLists.txt						\
 eigen-3.2.1/Eigen/Cholesky							\
 eigen-3.2.1/Eigen/QR								\
 eigen-3.2.1/Eigen/SuperLUSupport						\
 eigen-3.2.1/Eigen/Eigen							\
 eigen-3.2.1/Eigen/CholmodSupport						\
 eigen-3.2.1/Eigen/src/CholmodSupport/CholmodSupport.h				\
 eigen-3.2.1/Eigen/src/CholmodSupport/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h		\
 eigen-3.2.1/Eigen/src/IterativeLinearSolvers/IterativeSolverBase.h		\
 eigen-3.2.1/Eigen/src/IterativeLinearSolvers/CMakeLists.txt			\
 eigen-3.2.1/Eigen/src/IterativeLinearSolvers/IncompleteLUT.h			\
 eigen-3.2.1/Eigen/src/IterativeLinearSolvers/ConjugateGradient.h		\
 eigen-3.2.1/Eigen/src/IterativeLinearSolvers/BiCGSTAB.h			\
 eigen-3.2.1/Eigen/src/QR/ColPivHouseholderQR_MKL.h				\
 eigen-3.2.1/Eigen/src/QR/ColPivHouseholderQR.h					\
 eigen-3.2.1/Eigen/src/QR/FullPivHouseholderQR.h				\
 eigen-3.2.1/Eigen/src/QR/HouseholderQR.h					\
 eigen-3.2.1/Eigen/src/QR/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/QR/HouseholderQR_MKL.h					\
 eigen-3.2.1/Eigen/src/SparseCore/SparseBlock.h					\
 eigen-3.2.1/Eigen/src/SparseCore/SparsePermutation.h				\
 eigen-3.2.1/Eigen/src/SparseCore/SparseDiagonalProduct.h			\
 eigen-3.2.1/Eigen/src/SparseCore/SparseDenseProduct.h				\
 eigen-3.2.1/Eigen/src/SparseCore/SparseMatrixBase.h				\
 eigen-3.2.1/Eigen/src/SparseCore/CompressedStorage.h				\
 eigen-3.2.1/Eigen/src/SparseCore/AmbiVector.h					\
 eigen-3.2.1/Eigen/src/SparseCore/SparseUtil.h					\
 eigen-3.2.1/Eigen/src/SparseCore/SparseMatrix.h				\
 eigen-3.2.1/Eigen/src/SparseCore/MappedSparseMatrix.h				\
 eigen-3.2.1/Eigen/src/SparseCore/SparseCwiseBinaryOp.h				\
 eigen-3.2.1/Eigen/src/SparseCore/SparseCwiseUnaryOp.h				\
 eigen-3.2.1/Eigen/src/SparseCore/SparseSelfAdjointView.h			\
 eigen-3.2.1/Eigen/src/SparseCore/SparseTranspose.h				\
 eigen-3.2.1/Eigen/src/SparseCore/SparseFuzzy.h					\
 eigen-3.2.1/Eigen/src/SparseCore/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/SparseCore/ConservativeSparseSparseProduct.h		\
 eigen-3.2.1/Eigen/src/SparseCore/SparseSparseProductWithPruning.h		\
 eigen-3.2.1/Eigen/src/SparseCore/SparseView.h					\
 eigen-3.2.1/Eigen/src/SparseCore/SparseProduct.h				\
 eigen-3.2.1/Eigen/src/SparseCore/TriangularSolver.h				\
 eigen-3.2.1/Eigen/src/SparseCore/SparseColEtree.h				\
 eigen-3.2.1/Eigen/src/SparseCore/SparseVector.h				\
 eigen-3.2.1/Eigen/src/SparseCore/SparseDot.h					\
 eigen-3.2.1/Eigen/src/SparseCore/SparseTriangularView.h			\
 eigen-3.2.1/Eigen/src/SparseCore/SparseRedux.h					\
 eigen-3.2.1/Eigen/src/Householder/HouseholderSequence.h			\
 eigen-3.2.1/Eigen/src/Householder/Householder.h				\
 eigen-3.2.1/Eigen/src/Householder/BlockHouseholder.h				\
 eigen-3.2.1/Eigen/src/Householder/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h			\
 eigen-3.2.1/Eigen/src/Eigenvalues/GeneralizedEigenSolver.h			\
 eigen-3.2.1/Eigen/src/Eigenvalues/RealSchur.h					\
 eigen-3.2.1/Eigen/src/Eigenvalues/Tridiagonalization.h				\
 eigen-3.2.1/Eigen/src/Eigenvalues/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/Eigenvalues/HessenbergDecomposition.h			\
 eigen-3.2.1/Eigen/src/Eigenvalues/RealQZ.h					\
 eigen-3.2.1/Eigen/src/Eigenvalues/SelfAdjointEigenSolver_MKL.h			\
 eigen-3.2.1/Eigen/src/Eigenvalues/EigenSolver.h				\
 eigen-3.2.1/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h		\
 eigen-3.2.1/Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h			\
 eigen-3.2.1/Eigen/src/Eigenvalues/RealSchur_MKL.h				\
 eigen-3.2.1/Eigen/src/Eigenvalues/ComplexSchur.h				\
 eigen-3.2.1/Eigen/src/Eigenvalues/ComplexSchur_MKL.h				\
 eigen-3.2.1/Eigen/src/Eigenvalues/ComplexEigenSolver.h				\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_copy_to_ucol.h				\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU.h					\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_panel_bmod.h				\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_panel_dfs.h				\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_heap_relax_snode.h			\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_column_dfs.h				\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_Structs.h				\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_pivotL.h				\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_relax_snode.h				\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_Memory.h				\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLUImpl.h					\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_pruneL.h				\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_gemm_kernel.h				\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_SupernodalMatrix.h			\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_kernel_bmod.h				\
 eigen-3.2.1/Eigen/src/SparseLU/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_column_bmod.h				\
 eigen-3.2.1/Eigen/src/SparseLU/SparseLU_Utils.h				\
 eigen-3.2.1/Eigen/src/SuperLUSupport/SuperLUSupport.h				\
 eigen-3.2.1/Eigen/src/SuperLUSupport/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/CMakeLists.txt						\
 eigen-3.2.1/Eigen/src/Geometry/Transform.h					\
 eigen-3.2.1/Eigen/src/Geometry/Homogeneous.h					\
 eigen-3.2.1/Eigen/src/Geometry/OrthoMethods.h					\
 eigen-3.2.1/Eigen/src/Geometry/Rotation2D.h					\
 eigen-3.2.1/Eigen/src/Geometry/EulerAngles.h					\
 eigen-3.2.1/Eigen/src/Geometry/RotationBase.h					\
 eigen-3.2.1/Eigen/src/Geometry/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/Geometry/arch/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/Geometry/arch/Geometry_SSE.h				\
 eigen-3.2.1/Eigen/src/Geometry/ParametrizedLine.h				\
 eigen-3.2.1/Eigen/src/Geometry/Scaling.h					\
 eigen-3.2.1/Eigen/src/Geometry/AlignedBox.h					\
 eigen-3.2.1/Eigen/src/Geometry/Translation.h					\
 eigen-3.2.1/Eigen/src/Geometry/AngleAxis.h					\
 eigen-3.2.1/Eigen/src/Geometry/Hyperplane.h					\
 eigen-3.2.1/Eigen/src/Geometry/Umeyama.h					\
 eigen-3.2.1/Eigen/src/Geometry/Quaternion.h					\
 eigen-3.2.1/Eigen/src/StlSupport/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/StlSupport/StdVector.h					\
 eigen-3.2.1/Eigen/src/StlSupport/details.h					\
 eigen-3.2.1/Eigen/src/StlSupport/StdList.h					\
 eigen-3.2.1/Eigen/src/StlSupport/StdDeque.h					\
 eigen-3.2.1/Eigen/src/SPQRSupport/SuiteSparseQRSupport.h			\
 eigen-3.2.1/Eigen/src/SPQRSupport/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/PardisoSupport/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/PardisoSupport/PardisoSupport.h				\
 eigen-3.2.1/Eigen/src/Cholesky/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/Cholesky/LLT_MKL.h					\
 eigen-3.2.1/Eigen/src/Cholesky/LLT.h						\
 eigen-3.2.1/Eigen/src/Cholesky/LDLT.h						\
 eigen-3.2.1/Eigen/src/OrderingMethods/Amd.h					\
 eigen-3.2.1/Eigen/src/OrderingMethods/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/OrderingMethods/Eigen_Colamd.h				\
 eigen-3.2.1/Eigen/src/OrderingMethods/Ordering.h				\
 eigen-3.2.1/Eigen/src/LU/PartialPivLU_MKL.h					\
 eigen-3.2.1/Eigen/src/LU/FullPivLU.h						\
 eigen-3.2.1/Eigen/src/LU/PartialPivLU.h					\
 eigen-3.2.1/Eigen/src/LU/Determinant.h						\
 eigen-3.2.1/Eigen/src/LU/Inverse.h						\
 eigen-3.2.1/Eigen/src/LU/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/LU/arch/Inverse_SSE.h					\
 eigen-3.2.1/Eigen/src/LU/arch/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/Eigen2Support/VectorBlock.h				\
 eigen-3.2.1/Eigen/src/Eigen2Support/Minor.h					\
 eigen-3.2.1/Eigen/src/Eigen2Support/Geometry/Translation.h			\
 eigen-3.2.1/Eigen/src/Eigen2Support/Geometry/All.h				\
 eigen-3.2.1/Eigen/src/Eigen2Support/Geometry/CMakeLists.txt			\
 eigen-3.2.1/Eigen/src/Eigen2Support/Geometry/RotationBase.h			\
 eigen-3.2.1/Eigen/src/Eigen2Support/Geometry/Quaternion.h			\
 eigen-3.2.1/Eigen/src/Eigen2Support/Geometry/AlignedBox.h			\
 eigen-3.2.1/Eigen/src/Eigen2Support/Geometry/Transform.h			\
 eigen-3.2.1/Eigen/src/Eigen2Support/Geometry/Hyperplane.h			\
 eigen-3.2.1/Eigen/src/Eigen2Support/Geometry/ParametrizedLine.h		\
 eigen-3.2.1/Eigen/src/Eigen2Support/Geometry/Scaling.h				\
 eigen-3.2.1/Eigen/src/Eigen2Support/Geometry/Rotation2D.h			\
 eigen-3.2.1/Eigen/src/Eigen2Support/Geometry/AngleAxis.h			\
 eigen-3.2.1/Eigen/src/Eigen2Support/LeastSquares.h				\
 eigen-3.2.1/Eigen/src/Eigen2Support/Meta.h					\
 eigen-3.2.1/Eigen/src/Eigen2Support/Macros.h					\
 eigen-3.2.1/Eigen/src/Eigen2Support/Block.h					\
 eigen-3.2.1/Eigen/src/Eigen2Support/Cwise.h					\
 eigen-3.2.1/Eigen/src/Eigen2Support/TriangularSolver.h				\
 eigen-3.2.1/Eigen/src/Eigen2Support/LU.h					\
 eigen-3.2.1/Eigen/src/Eigen2Support/Lazy.h					\
 eigen-3.2.1/Eigen/src/Eigen2Support/CwiseOperators.h				\
 eigen-3.2.1/Eigen/src/Eigen2Support/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/Eigen2Support/Memory.h					\
 eigen-3.2.1/Eigen/src/Eigen2Support/SVD.h					\
 eigen-3.2.1/Eigen/src/Eigen2Support/MathFunctions.h				\
 eigen-3.2.1/Eigen/src/Eigen2Support/QR.h					\
 eigen-3.2.1/Eigen/src/SparseQR/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/SparseQR/SparseQR.h					\
 eigen-3.2.1/Eigen/src/PaStiXSupport/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/PaStiXSupport/PaStiXSupport.h				\
 eigen-3.2.1/Eigen/src/MetisSupport/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/MetisSupport/MetisSupport.h				\
 eigen-3.2.1/Eigen/src/Jacobi/Jacobi.h						\
 eigen-3.2.1/Eigen/src/Jacobi/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/UmfPackSupport/UmfPackSupport.h				\
 eigen-3.2.1/Eigen/src/UmfPackSupport/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/misc/SparseSolve.h					\
 eigen-3.2.1/Eigen/src/misc/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/misc/Image.h						\
 eigen-3.2.1/Eigen/src/misc/Kernel.h						\
 eigen-3.2.1/Eigen/src/misc/Solve.h						\
 eigen-3.2.1/Eigen/src/misc/blas.h						\
 eigen-3.2.1/Eigen/src/SparseCholesky/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/SparseCholesky/SimplicialCholesky.h			\
 eigen-3.2.1/Eigen/src/SparseCholesky/SimplicialCholesky_impl.h			\
 eigen-3.2.1/Eigen/src/plugins/CommonCwiseUnaryOps.h				\
 eigen-3.2.1/Eigen/src/plugins/MatrixCwiseBinaryOps.h				\
 eigen-3.2.1/Eigen/src/plugins/ArrayCwiseBinaryOps.h				\
 eigen-3.2.1/Eigen/src/plugins/CommonCwiseBinaryOps.h				\
 eigen-3.2.1/Eigen/src/plugins/BlockMethods.h					\
 eigen-3.2.1/Eigen/src/plugins/MatrixCwiseUnaryOps.h				\
 eigen-3.2.1/Eigen/src/plugins/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/plugins/ArrayCwiseUnaryOps.h				\
 eigen-3.2.1/Eigen/src/SVD/JacobiSVD_MKL.h					\
 eigen-3.2.1/Eigen/src/SVD/UpperBidiagonalization.h				\
 eigen-3.2.1/Eigen/src/SVD/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/SVD/JacobiSVD.h						\
 eigen-3.2.1/Eigen/src/Core/Ref.h						\
 eigen-3.2.1/Eigen/src/Core/DiagonalProduct.h					\
 eigen-3.2.1/Eigen/src/Core/VectorBlock.h					\
 eigen-3.2.1/Eigen/src/Core/Reverse.h						\
 eigen-3.2.1/Eigen/src/Core/Array.h						\
 eigen-3.2.1/Eigen/src/Core/Matrix.h						\
 eigen-3.2.1/Eigen/src/Core/MatrixBase.h					\
 eigen-3.2.1/Eigen/src/Core/util/ForwardDeclarations.h				\
 eigen-3.2.1/Eigen/src/Core/util/Meta.h						\
 eigen-3.2.1/Eigen/src/Core/util/Constants.h					\
 eigen-3.2.1/Eigen/src/Core/util/StaticAssert.h					\
 eigen-3.2.1/Eigen/src/Core/util/Macros.h					\
 eigen-3.2.1/Eigen/src/Core/util/BlasUtil.h					\
 eigen-3.2.1/Eigen/src/Core/util/XprHelper.h					\
 eigen-3.2.1/Eigen/src/Core/util/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/Core/util/Memory.h					\
 eigen-3.2.1/Eigen/src/Core/util/NonMPL2.h					\
 eigen-3.2.1/Eigen/src/Core/util/DisableStupidWarnings.h			\
 eigen-3.2.1/Eigen/src/Core/util/ReenableStupidWarnings.h			\
 eigen-3.2.1/Eigen/src/Core/util/MKL_support.h					\
 eigen-3.2.1/Eigen/src/Core/GenericPacketMath.h					\
 eigen-3.2.1/Eigen/src/Core/Block.h						\
 eigen-3.2.1/Eigen/src/Core/GlobalFunctions.h					\
 eigen-3.2.1/Eigen/src/Core/NumTraits.h						\
 eigen-3.2.1/Eigen/src/Core/NoAlias.h						\
 eigen-3.2.1/Eigen/src/Core/ProductBase.h					\
 eigen-3.2.1/Eigen/src/Core/Functors.h						\
 eigen-3.2.1/Eigen/src/Core/PlainObjectBase.h					\
 eigen-3.2.1/Eigen/src/Core/CoreIterators.h					\
 eigen-3.2.1/Eigen/src/Core/products/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/Core/products/CoeffBasedProduct.h			\
 eigen-3.2.1/Eigen/src/Core/products/SelfadjointProduct.h			\
 eigen-3.2.1/Eigen/src/Core/products/TriangularMatrixMatrix_MKL.h		\
 eigen-3.2.1/Eigen/src/Core/products/SelfadjointRank2Update.h			\
 eigen-3.2.1/Eigen/src/Core/products/SelfadjointMatrixMatrix.h			\
 eigen-3.2.1/Eigen/src/Core/products/GeneralMatrixMatrix.h			\
 eigen-3.2.1/Eigen/src/Core/products/GeneralMatrixMatrixTriangular_MKL.h	\
 eigen-3.2.1/Eigen/src/Core/products/TriangularMatrixVector.h			\
 eigen-3.2.1/Eigen/src/Core/products/TriangularSolverVector.h			\
 eigen-3.2.1/Eigen/src/Core/products/GeneralMatrixVector_MKL.h			\
 eigen-3.2.1/Eigen/src/Core/products/SelfadjointMatrixMatrix_MKL.h		\
 eigen-3.2.1/Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h		\
 eigen-3.2.1/Eigen/src/Core/products/Parallelizer.h				\
 eigen-3.2.1/Eigen/src/Core/products/GeneralMatrixMatrix_MKL.h			\
 eigen-3.2.1/Eigen/src/Core/products/SelfadjointMatrixVector_MKL.h		\
 eigen-3.2.1/Eigen/src/Core/products/GeneralBlockPanelKernel.h			\
 eigen-3.2.1/Eigen/src/Core/products/TriangularMatrixMatrix.h			\
 eigen-3.2.1/Eigen/src/Core/products/TriangularSolverMatrix.h			\
 eigen-3.2.1/Eigen/src/Core/products/TriangularMatrixVector_MKL.h		\
 eigen-3.2.1/Eigen/src/Core/products/SelfadjointMatrixVector.h			\
 eigen-3.2.1/Eigen/src/Core/products/GeneralMatrixVector.h			\
 eigen-3.2.1/Eigen/src/Core/products/TriangularSolverMatrix_MKL.h		\
 eigen-3.2.1/Eigen/src/Core/GeneralProduct.h					\
 eigen-3.2.1/Eigen/src/Core/CwiseNullaryOp.h					\
 eigen-3.2.1/Eigen/src/Core/Select.h						\
 eigen-3.2.1/Eigen/src/Core/CwiseBinaryOp.h					\
 eigen-3.2.1/Eigen/src/Core/SelfCwiseBinaryOp.h					\
 eigen-3.2.1/Eigen/src/Core/ArrayWrapper.h					\
 eigen-3.2.1/Eigen/src/Core/Stride.h						\
 eigen-3.2.1/Eigen/src/Core/CwiseUnaryOp.h					\
 eigen-3.2.1/Eigen/src/Core/IO.h						\
 eigen-3.2.1/Eigen/src/Core/EigenBase.h						\
 eigen-3.2.1/Eigen/src/Core/Assign.h						\
 eigen-3.2.1/Eigen/src/Core/MathFunctions.h					\
 eigen-3.2.1/Eigen/src/Core/Transpose.h						\
 eigen-3.2.1/Eigen/src/Core/ForceAlignedAccess.h				\
 eigen-3.2.1/Eigen/src/Core/arch/SSE/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/Core/arch/SSE/PacketMath.h				\
 eigen-3.2.1/Eigen/src/Core/arch/SSE/MathFunctions.h				\
 eigen-3.2.1/Eigen/src/Core/arch/SSE/Complex.h					\
 eigen-3.2.1/Eigen/src/Core/arch/NEON/PacketMath.h				\
 eigen-3.2.1/Eigen/src/Core/arch/NEON/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/Core/arch/NEON/Complex.h					\
 eigen-3.2.1/Eigen/src/Core/arch/Default/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/Core/arch/Default/Settings.h				\
 eigen-3.2.1/Eigen/src/Core/arch/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/Core/arch/AltiVec/Complex.h				\
 eigen-3.2.1/Eigen/src/Core/arch/AltiVec/PacketMath.h				\
 eigen-3.2.1/Eigen/src/Core/arch/AltiVec/CMakeLists.txt				\
 eigen-3.2.1/Eigen/src/Core/SelfAdjointView.h					\
 eigen-3.2.1/Eigen/src/Core/Swap.h						\
 eigen-3.2.1/Eigen/src/Core/MapBase.h						\
 eigen-3.2.1/Eigen/src/Core/Dot.h						\
 eigen-3.2.1/Eigen/src/Core/Fuzzy.h						\
 eigen-3.2.1/Eigen/src/Core/DenseCoeffsBase.h					\
 eigen-3.2.1/Eigen/src/Core/TriangularMatrix.h					\
 eigen-3.2.1/Eigen/src/Core/Replicate.h						\
 eigen-3.2.1/Eigen/src/Core/StableNorm.h					\
 eigen-3.2.1/Eigen/src/Core/CMakeLists.txt					\
 eigen-3.2.1/Eigen/src/Core/Random.h						\
 eigen-3.2.1/Eigen/src/Core/ReturnByValue.h					\
 eigen-3.2.1/Eigen/src/Core/BooleanRedux.h					\
 eigen-3.2.1/Eigen/src/Core/ArrayBase.h						\
 eigen-3.2.1/Eigen/src/Core/SolveTriangular.h					\
 eigen-3.2.1/Eigen/src/Core/Visitor.h						\
 eigen-3.2.1/Eigen/src/Core/VectorwiseOp.h					\
 eigen-3.2.1/Eigen/src/Core/BandMatrix.h					\
 eigen-3.2.1/Eigen/src/Core/DenseStorage.h					\
 eigen-3.2.1/Eigen/src/Core/NestByValue.h					\
 eigen-3.2.1/Eigen/src/Core/Assign_MKL.h					\
 eigen-3.2.1/Eigen/src/Core/Transpositions.h					\
 eigen-3.2.1/Eigen/src/Core/CwiseUnaryView.h					\
 eigen-3.2.1/Eigen/src/Core/CommaInitializer.h					\
 eigen-3.2.1/Eigen/src/Core/Diagonal.h						\
 eigen-3.2.1/Eigen/src/Core/DiagonalMatrix.h					\
 eigen-3.2.1/Eigen/src/Core/Map.h						\
 eigen-3.2.1/Eigen/src/Core/Flagged.h						\
 eigen-3.2.1/Eigen/src/Core/PermutationMatrix.h					\
 eigen-3.2.1/Eigen/src/Core/DenseBase.h						\
 eigen-3.2.1/Eigen/src/Core/Redux.h						\
 eigen-3.2.1/Eigen/LeastSquares							\
 eigen-3.2.1/COPYING.BSD							\
 eigen-3.2.1/COPYING.LGPL							\
 eigen-3.2.1/COPYING.MPL2							\
 eigen-3.2.1/COPYING.MINPACK							\
 eigen-3.2.1/INSTALL								\
 eigen-3.2.1/COPYING.README							\
 eigen-3.2.1/COPYING.GPL

if EXTERNAL_EIGEN
AM_CPPFLAGS = -I${EIGEN_PATH}
else
AM_CPPFLAGS = -I${top_srcdir}/${EIGEN_PATH}
endif

bin_PROGRAMS =
if BUILD_palinear
bin_PROGRAMS += palinear
endif

if BUILD_palogist
bin_PROGRAMS += palogist
endif

if BUILD_pacoxph
## Extra files for pacoxph
COXBASE = coxfit2 chinv2 cholesky2 chsolve2 dmatrix
COXSRC = $(COXBASE:=.c)
bin_PROGRAMS += pacoxph
endif

if BUILD_extractsnp
bin_PROGRAMS += extract-snp
endif


palinear_SOURCES = $(REGFILES) $(FVSRC) $(FVHEADERS) $(EIGENFILES)
palinear_CPPFLAGS = -DLINEAR $(AM_CPPFLAGS)

palogist_SOURCES = $(REGFILES) $(FVSRC) $(FVHEADERS) $(EIGENFILES)
palogist_CPPFLAGS = -DLOGISTIC $(AM_CPPFLAGS)

pacoxph_SOURCES = $(COXSRC) $(REGFILES) $(FVSRC) $(FVHEADERS)	\
 $(RHEADERS) survS.h survproto.h coxph_data.h coxph_data.cpp	\
 $(EIGENFILES)
pacoxph_CXXFLAGS = -I $(top_srcdir)/src/include $(AM_CXXFLAGS)
pacoxph_CPPFLAGS = -DCOXPH $(AM_CPPFLAGS)
pacoxph_CFLAGS = -I $(top_srcdir)/src/include $(AM_CFLAGS)

extract_snp_SOURCES = extract-snp.cpp $(FVSRC) $(FVHEADERS)

## Install these scripts in the bin directory as well:
dist_bin_SCRIPTS = probabel

## Install this R script in the examples directory
scriptdir = $(pkgdatadir)/scripts
dist_script_DATA = prepare_data.R extIDS.pl

# Copy the Eigen tar-ball to the ProbABEL tar-ball
EXTRA_DIST = 3.2.1.tar.bz2

## Install the config file
dist_sysconf_DATA = probabel_config.cfg.example

## Insert the correct values of several variables in the installed
## probabel script
install-exec-hook:
	$(SED) -i -e "s;probabel_config.cfg;@sysconfdir@/probabel_config.cfg;" \
		$(DESTDIR)$(bindir)/probabel
	$(SED) -i -e "s;\./;@bindir@/;" $(DESTDIR)$(bindir)/probabel
	$(SED) -i -e "s;PROBABEL_VERSION;@VERSION@;" $(DESTDIR)$(bindir)/probabel
	ln -s $(DESTDIR)$(bindir)/probabel $(DESTDIR)$(bindir)/probabel.pl

## Make the script executable
install-data-hook:
	chmod a+rx $(DESTDIR)/$(scriptdir)/extIDS.pl
