{-|
Module      : HABQTlib

This module contains functions for performing and simulating HABQT in Haskell.

Note: functions in this module simply call API from "HABQTlib.UnsafeAPI" after validating inputs.
-}
module HABQTlib
  ( TomState
  , TomFun
  , tomographyFun
  , simulatedTomography
  , streamResults
  ) where

import Control.Monad.State.Lazy
import Data.Validation
import HABQTlib.Data
import HABQTlib.MeasurementProcessing
import HABQTlib.UnsafeAPI
import Streaming (Of, Stream)
import qualified System.Random.MWC as MWC

validQMO ::
     QBitNum
  -> MHMCiter
  -> OptIter
  -> ( Validation [String] QBitNum
     , Validation [String] MHMCiter
     , Validation [String] OptIter)
validQMO nq' mi' oi' =
  let nq = validQBitN nq'
      mi = validMHMCiter mi'
      oi = validOptIter oi'
   in (nq, mi, oi)

-- | Given parameters such as output verbosity level and number of quantum
-- bits, set up the tomography function.
tomographyFun ::
     QBitNum -- ^ Number of quantum bits under tomography
  -> MHMCiter -- ^ Number of MHMC iterations to perform when resampling
  -> OptIter -- ^ Number of POVM optimisation steps to perform
  -> OutputVerb -- ^ Verbosity of stdout output
  -> MWC.GenIO -- ^ IO generator for variates from "System.Random.MWC"
  -> Validation [String] TomFun
tomographyFun nq' mi' oi' outv gen =
  let (nq, mi, oi) = validQMO nq' mi' oi'
   in tomographyFun' <$> nq <*> mi <*> oi <*> Success outv <*> Success gen

-- | Given a true state's density matrix and parameters, set up a simulation of
-- quantum tomography that outputs infidelity between mean estimates and true
-- state.
simulatedTomography ::
     DensityMatrix -- ^ True state's density matrix
  -> QBitNum -- ^ Number of quantum bits under tomography
  -> MHMCiter -- ^ Number of MHMC iterations to perform when resampling
  -> OptIter -- ^ Number of POVM optimisation steps to perform
  -> OutputVerb -- ^ Verbosity of stdout output
  -> MWC.GenIO -- ^ IO generator for variates from "System.Random.MWC"
  -> Validation [String] (StateT PurePOVM TomState Double)
simulatedTomography trueDM nq' mi' oi' outv gen =
  let dm = validDM trueDM
      (nq, mi, oi) = validQMO nq' mi' oi'
   in simulatedTomography' <$> dm <*> nq <*> mi <*> oi <*> Success outv <*>
      Success gen

-- | Stream simulated tomography results.
streamResults ::
     QBitNum -- ^ Number of quantum bits under tomography
  -> Rank -- ^ Rank of true state
  -> NumberOfParticles -- ^ Number of particles (per rank) to use for tomography
  -> MHMCiter -- ^ Number of MHMC iterations to perform when resampling
  -> OptIter -- ^ Number of POVM optimisation steps to perform
  -> OutputVerb -- ^ Verbosity of stdout output
  -> Validation [String] (Stream (Of Double) IO ())
streamResults nq' rank' pn' mi' oi' outv =
  let (nq, mi, oi) = validQMO nq' mi' oi'
      rank = validRank rank'
      pn = validPartNum pn'
   in streamResults' <$> nq <*> rank <*> pn <*> mi <*> oi <*> Success outv
