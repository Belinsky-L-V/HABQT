{-|
Module      : HABQTlib.UnsafeAPI

This module contains functions for performing and simulating HABQT in Haskell.

__Caution__: functions in this module perform no input validation and are partial. For a safe API refer to "HABQTlib".
-}
module HABQTlib.UnsafeAPI where

import Control.Applicative (liftA2)
import Control.Monad (replicateM)
import Control.Monad.State.Lazy
import Data.Maybe (fromJust)
import qualified Data.Vector as V
import HABQTlib.Data
import HABQTlib.Data.Particle
import HABQTlib.MeasurementProcessing
import HABQTlib.RandomStates
import qualified Numeric.LinearAlgebra as LA
import Streaming
import qualified Streaming.Prelude as S
import qualified System.Random.MWC as MWC
import Text.Printf (printf)

-- | Tomography keeps track of the particle hierarchy and list of previous
-- measurement results, IO is used for verbose output and assorted random state
-- generation.
type TomState = StateT (ParticleHierarchy, [PureStateVector]) IO

-- | Tomography function takes a measurement result and returns state-dependent
-- Bayesian mean estimate of state and the optimal next POVM to perform.
type TomFun = PureStateVector -> TomState (DensityMatrix, PurePOVM)

-- | Given parameters such as output verbosity level and number of quantum
-- bits, set up the tomography function.
tomographyFun' ::
     QBitNum -- ^ Number of quantum bits under tomography
  -> MHMCiter -- ^ Number of MHMC iterations to perform when resampling
  -> OptIter -- ^ Number of POVM optimisation steps to perform
  -> OutputVerb -- ^ Verbosity of stdout output
  -> MWC.GenIO -- ^ IO generator for variates from "System.Random.MWC"
  -> TomFun
tomographyFun' nq mi oi outv gen nextResult = do
  (ph, ms) <- get
  let nextPH = updateParticleHierarchy nextResult ph
      dim = LA.rows . getStateVector $ nextResult
      effectiveSizes =
        V.map (liftA2 (/) effectiveSize (fromIntegral . ptsNumber)) nextPH
      ra = ResampleArgs outv gen dim mi
      resampleC es pts =
        if es < 0.5
          then resample ra (nextResult : ms) pts
          else return pts
  nextPH' <- lift $ V.zipWithM resampleC effectiveSizes nextPH
  sv0s <- liftIO $ replicateM nq (genPureSV 2)
  let nextEstimate = getMixedEstimate nextPH'
      nextPOVM = optimiseSingleQbPOVM oi sv0s nextPH'
  put (nextPH', nextResult : ms)
  return (nextEstimate, nextPOVM)

-- | Given a true state's density matrix and parameters, set up a simulation of
-- quantum tomography that outputs infidelity between mean estimates and true
-- state.
simulatedTomography' ::
     DensityMatrix -- ^ True state's density matrix
  -> QBitNum -- ^ Number of quantum bits under tomography
  -> MHMCiter -- ^ Number of MHMC iterations to perform when resampling
  -> OptIter -- ^ Number of POVM optimisation steps to perform
  -> OutputVerb -- ^ Verbosity of stdout output
  -> MWC.GenIO -- ^ IO generator for variates from "System.Random.MWC"
  -> StateT PurePOVM TomState Double
simulatedTomography' trueDM nq mi oi outv gen = do
  povm <- get
  nextResult <- liftIO $ simulateMeasuremet trueDM povm gen
  (nextEstimate, nextPOVM) <- lift $ tomographyFun' nq mi oi outv gen nextResult
  (nextPH, _) <- lift get
  put nextPOVM
  let fid = fidelityDM trueDM nextEstimate
  when (outv > NoOutput) . liftIO $ do
    let dim = 2 ^ nq
        rankFids =
          V.map
            (fidelityDM trueDM . snd . getWDM . reduceParticlesToMean)
            nextPH
        weightsAndFids =
          V.zip3 (V.enumFromN (1 :: Rank) dim) (V.map ptsWeight nextPH) rankFids
    putStrLn ""
    V.mapM_
      (\(a, b, c) ->
         printf "(Rank: %4d, Weight: %10.9f, Fidelity: %10.9f)\n" a b c)
      weightsAndFids
  return $ 1 - fid

-- | Stream simulated tomography results.
streamResults' ::
     QBitNum -- ^ Number of quantum bits under tomography
  -> Rank -- ^ Rank of true state
  -> NumberOfParticles -- ^ Number of particles (per rank) to use for tomography
  -> MHMCiter -- ^ Number of MHMC iterations to perform when resampling
  -> OptIter -- ^ Number of POVM optimisation steps to perform
  -> OutputVerb -- ^ Verbosity of stdout output
  -> Stream (Of Double) IO ()
streamResults' nq rank pn mi oi outv = do
  let dim = 2 ^ nq
  trueDM <- liftIO $ genDM dim rank
  ph <- liftIO $ initialiseParticleHierarchy dim pn
  gen <- liftIO MWC.createSystemRandom
  rPOVM <-
    liftIO $
    productPOVM <$>
    replicateM nq (mkAntipodalPOVM . fromJust . svToAngles <$> genPureSV 2)
  let tomS = S.repeatM (simulatedTomography' trueDM nq mi oi outv gen)
      tomS' = evalStateT (distribute tomS) rPOVM
      initInfid = 1 - fidelityDM trueDM (getMixedEstimate ph)
  S.yield initInfid
  evalStateT (distribute tomS') (ph, [])
