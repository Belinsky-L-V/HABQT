{-# LANGUAGE RecordWildCards #-}

{-|
Module      :  HABQTlib.Data.Particle

Data structures and functions that deal with storing and processing particle
hierarchies.

/Warning/: functions in this module assume that the 'ptsParticles' is non-empty
and 'NumberOfParticles', 'Dim', and 'Rank' are positive, no validation is
performed. If you use them directly, instead of employing API from
"HABQTlib", you must ensure those assumptions hold.
-}
module HABQTlib.Data.Particle
  ( Particles(..)
  , genParticles
  , updateParticles
  , ParticleHierarchy
  , initialiseParticleHierarchy
  , updateParticleHierarchy
  , getMixedEstimate
  , foldOverPts
  , reduceParticlesToMean
  , effectiveSize
  , ResampleArgs(..)
  , resampleMultinom
  , resample
  , ecdf
  , icdf
  , nudgeParticle
  ) where

import Control.Monad (when)
import Control.Newtype.Generics (over)
import Data.Bool.HT (if', select)
import qualified Data.Vector as V
import HABQTlib.Data
import HABQTlib.RandomStates
import Numeric.LinearAlgebra (Complex(..))
import qualified Numeric.LinearAlgebra as LA
import qualified System.Random.MWC as MWC
import Text.Printf (printf)

-- | A vector of weighed density matrices is stored along with their rank and
-- number. 'ptsWeight' corresponds to the collective weight of particles of
-- rank 'ptsRank' in the hierarchical model, it is not the sum of individual
-- weights of particles (that is normalised to unity after every update).
data Particles = Particles
  { ptsRank :: Rank
  , ptsWeight :: Weight
  , ptsNumber :: NumberOfParticles
  , ptsParticles :: V.Vector WeighedDensityMatrix
  } deriving (Show)

-- | Particle hierarchy is described by a vector of 'Particles'.
type ParticleHierarchy = V.Vector Particles

-- | Generates random particles (according to induced measure).
genParticles :: Dim -> Rank -> NumberOfParticles -> IO Particles
genParticles d r n =
  let w = 1 / fromIntegral n
   in Particles r 1 n . fmap (mkWDM w) <$> V.replicateM n (genDM d r)

-- | Summarise particles to a mean estimate (weighed by the corresponding
-- hierarchical weight of the rank).
reduceParticlesToMean :: Particles -> WeighedDensityMatrix
reduceParticlesToMean Particles {..} =
  let wdm = V.foldl1' (<+>) ptsParticles
      wdmw = over WeighedDensityMatrix (\(w, dm) -> (w * ptsWeight, dm)) wdm
   in truncateRank ptsRank wdmw

-- | Map density matrices, combine them with their weights, and then perform a
-- (strict left) fold of results.
foldOverPts ::
     (DensityMatrix -> a) -- ^ function to map over density matrices
  -> (Weight -> a -> b) -- ^ function to combine weights with results of mapping
  -> (c -> b -> c) -- ^ fold funciton
  -> c -- ^ seed value for folding
  -> Particles
  -> c
foldOverPts f wf fld z Particles {..} =
  let wm (WeighedDensityMatrix (w, dm)) = wf w (f dm)
   in V.foldl' (\l r -> fld l (wm r)) z ptsParticles

fullDataLogLikelihood :: [PureStateVector] -> DensityMatrix -> Double
fullDataLogLikelihood vs dm =
  let lps = map (log . (`pureStateLikelihood` dm)) vs
   in sum lps

-- | Given a measurement result, perform a Bayesian update over the particles.
updateParticles :: PureStateVector -> Particles -> Particles
updateParticles sv pts@Particles {..} =
  let updateF :: WeighedDensityMatrix -> WeighedDensityMatrix
      updateF (WeighedDensityMatrix (w, dm)) = WeighedDensityMatrix (wnew, dm)
        where
          wnew = w * pureStateLikelihood sv dm
      upts = V.map updateF ptsParticles
      uw = V.foldl' (\acc (WeighedDensityMatrix (w, _)) -> acc + w) 0 upts
      npts = V.map (over WeighedDensityMatrix (\(w, dm) -> (w / uw, dm))) upts
   in pts {ptsWeight = ptsWeight * uw, ptsParticles = npts}

-- | Helper function that generates random particles of each applicable rank.
initialiseParticleHierarchy :: Dim -> NumberOfParticles -> IO ParticleHierarchy
initialiseParticleHierarchy d n = V.generateM d (\r -> genParticles d (r + 1) n)

-- | Given a measurement result, update all particles, then normalise resulting
-- hierarchical weights to sum to unity.
updateParticleHierarchy ::
     PureStateVector -> ParticleHierarchy -> ParticleHierarchy
updateParticleHierarchy sv ph =
  let uph = V.map (updateParticles sv) ph
      wgts = V.map ptsWeight uph
      nwgts = V.map (/ V.sum wgts) wgts
   in V.zipWith (\x w -> x {ptsWeight = w}) uph nwgts

-- | Summarise the whole particle hierarchy into one mean Bayesian estimate.
getMixedEstimate :: ParticleHierarchy -> DensityMatrix
getMixedEstimate ph =
  let rankEstimates = V.map reduceParticlesToMean ph
      WeighedDensityMatrix (_, result) = V.foldl1' (<+>) rankEstimates
   in result

-- | Calculate the effective sample size of particles (weights don’t
-- necessarily have to be normalised).
effectiveSize :: Particles -> Double
effectiveSize Particles {..} =
  let ss = V.sum . V.map ((^ (2 :: Int)) . fst . getWDM) $ ptsParticles
      wa = V.foldl' (flip ((+) . fst . getWDM)) 0 ptsParticles
   in wa ^ (2 :: Int) / ss

-- | Nudges a particle by mixing the state together with some randomly
-- generated pure state. Relative weight of the random component determines how
-- “close” a nudged particle is to the original one.
nudgeParticle ::
     Dim
  -> Weight -- ^ Relative weight (from 0 to 1) of random component
  -> WeighedDensityMatrix
  -> IO WeighedDensityMatrix
nudgeParticle dim weightFraction (WeighedDensityMatrix (w, dm)) = do
  DensityMatrix nudgeDM <- svToDM <$> genPureSV dim
  let dmw = LA.scale (1 - (weightFraction :+ 0)) (getDensityMatrix dm)
      dmwn = LA.scale (weightFraction :+ 0) nudgeDM
  return $ WeighedDensityMatrix (w, DensityMatrix $ dmw + dmwn)

-- | Calculates values of empirical distribution function at data points.
ecdf :: V.Vector WeighedDensityMatrix -> V.Vector Double
ecdf = V.postscanl' (+) 0 . V.map (fst . getWDM)

-- | /O(log n)/ Given a non-empty sorted vector (typically an empirical cdf
-- evaluated at data points returned by ecdf) and a real number return the
-- (0-based) index of the least element of vector which is greater or equal to
-- the given real number (or the index of the last element, in case there is no
-- element smaller than the argument).
icdf :: V.Vector Double -> Double -> Int
icdf cdf x =
  let tIdx = V.length cdf - 1
      go (lIdx, hIdx) =
        let mIdx =
              truncate $ ((fromIntegral lIdx :: Double) + fromIntegral hIdx) / 2
         in select
              (go (lIdx, mIdx))
              [ (lIdx == hIdx, lIdx)
              , (lIdx + 1 == hIdx, if' (cdf V.! lIdx > x) lIdx hIdx)
              , (cdf V.! mIdx < x, go (mIdx, hIdx))
              ]
   in select (go (0, tIdx)) [(x <= V.head cdf, 0), (x > V.last cdf, tIdx)]

-- | Multinomial resampling of particle vector, which equalises weights of
-- particles.
resampleMultinom :: MWC.GenIO -> Particles -> IO Particles
resampleMultinom gen pts@Particles {..} = do
  us <- MWC.uniformVector gen ptsNumber
  let cdf = ecdf ptsParticles
      idxs = V.map (icdf cdf) us
      w = 1 / fromIntegral ptsNumber
      pointR = over WeighedDensityMatrix (\(_, dm) -> (w, dm))
      resampled = V.map (ptsParticles V.!) idxs
      normed = V.map pointR resampled
  return pts {ptsParticles = normed}

mhmcStep ::
     MWC.GenIO
  -> Dim
  -> Double
  -> [PureStateVector]
  -> Particles
  -> IO (Double, Particles)
mhmcStep gen dim rw ms pts@Particles {..} = do
  let cr wdm wdm' =
        exp
          (fullDataLogLikelihood ms (snd . getWDM $ wdm') -
           fullDataLogLikelihood ms (snd . getWDM $ wdm))
  newParticles <-
    V.mapM (fmap (truncateRank ptsRank) . nudgeParticle dim rw) ptsParticles
  us <- V.replicateM ptsNumber (MWC.uniform gen :: IO Double)
  let crs = V.zipWith cr ptsParticles newParticles
      change = V.zipWith (<=) us crs
      accRate =
        (fromIntegral . V.length . V.filter id) change / fromIntegral ptsNumber
      rwdms = V.zipWith3 if' change newParticles ptsParticles
      final = pts {ptsParticles = rwdms}
  return (accRate, final)

resampleMHMC ::
     ResampleArgs
  -> DensityMatrix
  -> Double
  -> Int
  -> [PureStateVector]
  -> Particles
  -> IO Particles
resampleMHMC ra@ResampleArgs {..} estimate wr iter mts pts = do
  (accRate, resampled) <- mhmcStep argGen argDim wr mts pts
  when (argOut == FullOutput) $
    printf
      "(Weight of new particle: %8.3g, MHMC acceptance rate: %8.3g)\n"
      wr
      accRate
  let (iter', wr') =
        select
          (iter + 1, wr)
          [ (accRate < 1e-2, (0, wr * 0.25))
          , (accRate < 1e-1, (0, wr * 0.5))
          , (iter < argMinIter, (iter + 1, wr))
          , (accRate < 0.33, (0, wr * 0.5))
          ]
  if iter > argMinIter
    then return resampled
    else resampleMHMC ra estimate wr' iter' mts resampled

-- | Arguments for the resampling function.
data ResampleArgs = ResampleArgs
  { argOut :: OutputVerb
  , argGen :: MWC.GenIO
  , argDim :: Dim
  , argMinIter :: MHMCiter
  }

-- | Resample particles. First does one multinomial step that equalises the
-- weights, then performs MHMC iterations adaptively refining the proposal
-- distribution based on acceptance rate. 'argMinIter' iterations are performed
-- for proposal distributions with significant acceptance rates.
resample :: ResampleArgs -> [PureStateVector] -> Particles -> IO Particles
resample ra@ResampleArgs {..} mts pts@Particles {..} = do
  let estimate = getMixedEstimate . V.singleton $ pts
      nudgeW = 0.95
  when (argOut == FullOutput) $ do
    putStrLn ""
    putStrLn $ "resampling rank " ++ show ptsRank
    putStrLn ""
  rm <- resampleMultinom argGen pts
  resampleMHMC ra estimate nudgeW 0 mts rm
