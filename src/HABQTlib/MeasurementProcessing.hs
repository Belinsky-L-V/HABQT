{-|
Module      : HABQTlib.MeasurementProcessing

Functions that deal with optimising and simulating measurements that take place
during tomography.
-}
module HABQTlib.MeasurementProcessing
  ( PurePOVM
  , SingleQbParam
  , measurementProbs
  , svToAngles
  , blochAnglesToSV
  , mkAntipodalPOVM
  , productPOVM
  , simulateMeasuremet
  , optimiseSingleQbPOVM
  ) where

import Control.Applicative (liftA2)
import Control.Newtype.Generics (unpack)
import Data.Complex (mkPolar, polar)
import Data.List (unfoldr)
import Data.Maybe (fromJust)
import qualified Data.Vector as V
import HABQTlib.Data
import HABQTlib.Data.Particle
import qualified Numeric.GSL as GSL
import Numeric.LinearAlgebra (Complex((:+)))
import qualified Numeric.LinearAlgebra as LA
import qualified System.Random.MWC as MWC
import System.Random.MWC.Distributions (categorical)

-- | A POVM consisting of projections onto pure states.
type PurePOVM = V.Vector WeighedPureStateVector

-- | Spherical coordinates of a pure single qubit state on Bloch sphere.
type SingleQbParam = (Double, Double)

-- | Probabilities to measure elements of POVM when performing the measurement
-- over a mixed state.
measurementProbs :: PurePOVM -> DensityMatrix -> V.Vector Double
measurementProbs povm dm = V.map (weighedProb dm) povm

weighedProb :: DensityMatrix -> WeighedPureStateVector -> Double
weighedProb dm (WeighedPureStateVector (w, sv)) = w * pureStateLikelihood sv dm

entropy :: V.Vector Double -> Double
entropy = negate . V.sum . V.map (\p -> p * log p)

povmPointEntropy :: PurePOVM -> DensityMatrix -> Double
povmPointEntropy povm dm = entropy (measurementProbs povm dm)

povmMeanEntropy :: PurePOVM -> ParticleHierarchy -> Double
povmMeanEntropy povm ph =
  let meanEnt = foldOverPts (povmPointEntropy povm) (*) (+) 0
      totalWeight = V.sum . V.map ptsWeight $ ph
      rawEntropy = V.sum . V.map (liftA2 (*) ptsWeight meanEnt) $ ph
   in rawEntropy / totalWeight

-- | Recover a state vector from spherical coordinates.
blochAnglesToSV :: SingleQbParam -> PureStateVector
blochAnglesToSV (th, phi) =
  let z = LA.fromList [1, 0]
      o = LA.fromList [0, 1]
   in PureStateVector . LA.fromColumns . pure $
      LA.scalar (cos (th / 2) :+ 0) * z +
      LA.scalar (mkPolar (sin (th / 2)) phi) * o

-- | Return spherical coordinates of a single qubit pure state (on Bloch
-- sphere).
svToAngles :: PureStateVector -> Maybe SingleQbParam
svToAngles (PureStateVector sv) =
  let s = LA.size sv
      (m0, ph0) = polar $ LA.atIndex sv (0, 0)
      (_, ph1) = polar $ LA.atIndex sv (1, 0)
      ph = ph1 - ph0
      th = 2 * acos m0
   in if s == (2, 1)
        then Just (th, ph)
        else Nothing

-- | Given spherical coordinates, construct a POVM from the given vector and
-- one orthogonal to it.
mkAntipodalPOVM :: SingleQbParam -> PurePOVM
mkAntipodalPOVM c@(th, phi) =
  let sv = blochAnglesToSV c
      sv' = blochAnglesToSV (pi - th, phi + pi)
   in V.fromList $ WeighedPureStateVector <$> [(1, sv), (1, sv')]

reshapeList :: Int -> [a] -> [[a]]
reshapeList n =
  unfoldr
    (\b ->
       if length b < n
         then Nothing
         else Just (splitAt n b))

listToPairs :: [a] -> [(a, a)]
listToPairs = fmap (\(a:[b]) -> (a, b)) . reshapeList 2

pairToList :: (a, a) -> [a]
pairToList (a, b) = [a, b]

-- | Given a list of POVM measurements on sub-systems, construct a POVM over
-- the composite system that includes all of them.
productPOVM :: [PurePOVM] -> PurePOVM
productPOVM sqbPovms =
  let pr ::
           WeighedPureStateVector
        -> WeighedPureStateVector
        -> WeighedPureStateVector
      pr wsv0 wsv1 =
        WeighedPureStateVector (w0 * w1, PureStateVector $ LA.kronecker sv0 sv1)
        where
          up = fmap unpack . unpack
          (w0, sv0) = up wsv0
          (w1, sv1) = up wsv1
   in V.foldl1' (liftA2 pr) $ V.fromList sqbPovms

-- | Approximate most informative separable POVM over a composite system of
-- quantum bits, given a list of single-qubit starting points and a particle
-- distribution.
optimiseSingleQbPOVM ::
     OptIter -- ^ Number of optimisation steps to perform
  -> [PureStateVector] -- ^ single qubit initial states
  -> ParticleHierarchy
  -> PurePOVM
optimiseSingleQbPOVM iter sv0s ph =
  let nq = length sv0s
      method = GSL.NMSimplex2
      precision = 1e-6
      iterations = iter
      initialBox = replicate (2 * nq) (pi / 2)
      estimateDM = getMixedEstimate ph
      obj params =
        negate $ povmPointEntropy povm estimateDM - povmMeanEntropy povm ph
        where
          povm = productPOVM . fmap mkAntipodalPOVM . listToPairs $ params
      start = concatMap (pairToList . fromJust . svToAngles) sv0s
      result = GSL.minimize method precision iterations initialBox obj start
   in productPOVM . fmap mkAntipodalPOVM . listToPairs . fst $ result

-- | Simulate a POVM over a mixed state and return the state vector on which
-- the projection was obtained.
simulateMeasuremet ::
     DensityMatrix -> PurePOVM -> MWC.GenIO -> IO PureStateVector
simulateMeasuremet dm povm gen = do
  let probs = measurementProbs povm dm
      svs = V.map (\(WeighedPureStateVector (_, sv)) -> sv) povm
  idx <- categorical probs gen
  return $ svs V.! idx
