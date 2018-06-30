{-# LANGUAGE DeriveGeneric #-}

{-|
Module      : HABQTlib.Data

This module contains data types and helper functions for working with quantum
state vectors and density matrices.
-}
module HABQTlib.Data
  ( Dim
  , Rank
  , NumberOfParticles
  , QBitNum
  , Weight
  , MHMCiter
  , OptIter
  , OutputVerb(..)
  , DensityMatrix(..)
  , truncateRank
  , PureStateVector(..)
  , pureStateLikelihood
  , svToDM
  , WeighedDensityMatrix(..)
  , mkWDM
  , mkWDM1
  , WeighedPureStateVector(..)
  , (<+>)
  , fidelity
  , fidelityDM
  , validSV
  , validDM
  , validPartNum
  , validRank
  , validMHMCiter
  , validOptIter
  , validQBitN
  ) where

import Control.Newtype.Generics (Newtype, unpack)
import Data.Bool.HT (select)
import Data.Complex
import Data.Validation
import GHC.Generics (Generic)
import qualified Numeric.LinearAlgebra as LA
import Numeric.LinearAlgebra (Matrix)

-- | Dimension of Hilbert space.
type Dim = Int

-- | Rank of mixed state.
type Rank = Int

-- | Number of particles per rank.
type NumberOfParticles = Int

-- | Number of quantum bits.
type QBitNum = Int

-- | Number of MHMC iterations to perform when resampling.
type MHMCiter = Int

-- | Number of optimisation steps to perform when searching for optimal
-- measurment.
type OptIter = Int

-- | Weight associated with a particle.
type Weight = Double

-- | Output verbosity settings.
data OutputVerb
  = NoOutput -- ^ No stdout output
  | FidOutput -- ^ Only output fidelities and weights of hierarchical mean estimates
  | FullOutput -- ^ Full output, including resampling diagnostic information
  deriving (Eq, Show, Ord)

-- | Density matrix are stored as hmatrix matrices of complex doubles.
newtype DensityMatrix = DensityMatrix
  { getDensityMatrix :: Matrix (Complex Double)
  } deriving (Eq, Show, Generic)

-- | Pure state vectors are stored as hmatrix matrices of complex doubles.
-- Such matrices only have one column.
newtype PureStateVector = PureStateVector
  { getStateVector :: Matrix (Complex Double)
  } deriving (Eq, Show, Generic)

instance Newtype DensityMatrix

instance Newtype PureStateVector

-- | Check whether a pure state vector is properly normed.
validSV :: PureStateVector -> Validation [String] PureStateVector
validSV =
  validate
    ["State vector must have unit norm."]
    (\x -> abs (1 - LA.norm_2 (unpack x)) < 1e-12)

-- | Verify that density matrix is Hermitian and has trace 1.
validDM :: DensityMatrix -> Validation [String] DensityMatrix
validDM =
  let traceU :: LA.Matrix (Complex Double) -> Bool
      traceU dm =
        (abs (1 - (magnitude . LA.sumElements . LA.takeDiag) dm) < 1e-12)
      hermU dm = (LA.norm_2 (dm - LA.tr dm) < 1e-6)
      both = ((&&) <$> traceU <*> hermU) . unpack
      dmM = ["Density matrix must be Hermitian and have trace of 1."]
   in validate dmM both

qnM :: [String]
qnM = ["Number of quantum bits must be a positive integer."]

-- | Verify that number of quantum bits is positive.
validQBitN :: QBitNum -> Validation [String] QBitNum
validQBitN = validate qnM (> 0)

miM :: [String]
miM = ["Number of MHMC iterations must be a positive integer."]

-- | Verify that number of MHMC iterations is a positive integer.
validMHMCiter :: MHMCiter -> Validation [String] MHMCiter
validMHMCiter = validate miM (> 0)

pnM :: [String]
pnM = ["Number of particles per rank must be a positive integer."]

-- | Verify that particle number is a positive integer.
validPartNum :: NumberOfParticles -> Validation [String] NumberOfParticles
validPartNum = validate pnM (> 0)

rM :: [String]
rM = ["Rank must be a positive integer."]

-- | Verify that rank is a positive integer. Setting rank to be higher than
-- the dimension of space creates poinless performance overhead, but isn't
-- prevented by validation.
validRank :: Rank -> Validation [String] Rank
validRank = validate rM (> 0)

oiM :: [String]
oiM = ["Number of POVM optimisation iterations must be a positive integer."]

-- | Verify that number of POVM optimisation iterations is positive.
validOptIter :: OptIter -> Validation [String] OptIter
validOptIter = validate oiM (> 0)

-- | Weighed density matrix where weight is stored separately as first
-- coordinate of a tuple.
newtype WeighedDensityMatrix = WeighedDensityMatrix
  { getWDM :: (Weight, DensityMatrix)
  } deriving (Eq, Show, Generic)

-- | A shorter alias for curried 'WeighedDensityMatrix' constructor.
mkWDM :: Weight -> DensityMatrix -> WeighedDensityMatrix
mkWDM w dm = WeighedDensityMatrix (w, dm)

-- | Alias for @'mkWDM' 1@.
mkWDM1 :: DensityMatrix -> WeighedDensityMatrix
mkWDM1 = mkWDM 1

-- | Weighed state vector where weight is stored separately as first coordinate
-- of a tuple.
newtype WeighedPureStateVector = WeighedPureStateVector
  { getWSV :: (Weight, PureStateVector)
  } deriving (Eq, Show, Generic)

instance Newtype WeighedDensityMatrix

instance Newtype WeighedPureStateVector

-- | Fidelity (probability of measurement) between pure states.
fidelity :: PureStateVector -> PureStateVector -> Double
fidelity (PureStateVector sv1) (PureStateVector sv2) =
  let ips = magnitude (LA.atIndex (LA.tr sv1 LA.<> sv2) (0, 0)) ^ (2 :: Int)
   in select ips [(ips < 0, 0), (ips > 1, 1)]

-- | Fidelity (probability of measurement) between mixed states.
fidelityDM :: DensityMatrix -> DensityMatrix -> Double
fidelityDM (DensityMatrix dm1) (DensityMatrix dm2) =
  let (u1, s1) = LA.leftSV dm1
      (u2, s2) = LA.leftSV dm2
      ss1 = LA.cmap (\x -> sqrt x :+ 0) s1
      ss2 = LA.cmap (\x -> sqrt x :+ 0) s2
   in LA.norm_nuclear
        (u1 LA.<> LA.diag ss1 LA.<> LA.tr u1 LA.<> u2 LA.<> LA.diag ss2 LA.<>
         LA.tr u2) ^
      (2 :: Int)

-- | Calculate the density matrix of a given pure state.
svToDM :: PureStateVector -> DensityMatrix
svToDM (PureStateVector sv) = DensityMatrix $ sv LA.<> LA.tr sv

infix 8 <+>

-- | Given two weighed density matrixes, compute their mixture. Associative
-- operation.
(<+>) :: WeighedDensityMatrix -> WeighedDensityMatrix -> WeighedDensityMatrix
wdm0 <+> wdm1 =
  WeighedDensityMatrix
    (w0 + w1, DensityMatrix $ LA.scale c0 dm0 + LA.scale c1 dm1)
  where
    up = fmap unpack . unpack
    (w0, dm0) = up wdm0
    (w1, dm1) = up wdm1
    cs = LA.fromList [w0 :+ 0, w1 :+ 0]
    csn = LA.scale (1 / LA.sumElements cs) cs
    c0 = csn LA.! 0
    c1 = csn LA.! 1

-- | Probability of obtaining a measurement result when projecting a system in
-- mixed state determined by a density matrix onto a pure state.
pureStateLikelihood :: PureStateVector -> DensityMatrix -> Double
pureStateLikelihood (PureStateVector sv) (DensityMatrix dm) =
  let p =
        LA.magnitude . LA.sumElements . LA.takeDiag $ LA.tr sv LA.<> dm LA.<> sv
   in select p [(p < 0, 0), (p > 1, 1)]

-- | Set smallest eigenvalues of a weighed density matrix to zero until
-- specified rank is reached.
truncateRank :: Rank -> WeighedDensityMatrix -> WeighedDensityMatrix
truncateRank targetRank (WeighedDensityMatrix (w, DensityMatrix dm)) =
  let (u, s, _) = LA.svd dm
      st = LA.real $ LA.subVector 0 targetRank s
      ut = u LA.¿ [0 .. (targetRank - 1)]
      stn = LA.scale (1 / LA.sumElements st) st
   in WeighedDensityMatrix
        (w, DensityMatrix $ ut LA.<> LA.diag stn LA.<> LA.tr ut)
