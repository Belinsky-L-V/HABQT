module StateGenTests
  ( testStateGen
  , testStateArb
  , arbPropPure
  , arbPropDensityMatrix
  ) where

import Data.Complex (magnitude)
import HABQTlib.Data
import HABQTlib.RandomStates
import Numeric.LinearAlgebra (norm_2, sumElements, takeDiag, tr)
import qualified System.Random.MWC as MWC
import qualified Test.QuickCheck as QC
import Test.QuickCheck (Property)
import TestHelpers

class StateProp a where
  traceHermProp :: MWC.GenIO -> a -> Property

dimRankIO :: MWC.GenIO -> Dim -> IO (Dim, Rank)
dimRankIO gen ub = do
  dim <- MWC.uniformR (1, ub) gen
  r <- MWC.uniformR (1, dim) gen
  return (dim, r)

instance StateProp DensityMatrix where
  traceHermProp gen _ =
    QC.ioProperty $ do
      (dim, r) <- dimRankIO gen 100
      dm <- getDensityMatrix <$> genDM dim r
      return $
        (abs (1 - (magnitude . sumElements . takeDiag) dm) < 1e-12) &&
        (norm_2 (dm - tr dm) < 1e-12)

instance StateProp PureStateVector where
  traceHermProp gen _ =
    QC.ioProperty $ do
      (dim, _) <- dimRankIO gen 100
      sv <- getStateVector <$> genPureSV dim
      return $ abs (1 - norm_2 sv) < 1e-12

testStateGen :: IO ()
testStateGen = do
  gen <- MWC.createSystemRandom
  putStrLn "Testing density matrix generation:"
  QC.quickCheckWith
    QC.stdArgs {QC.maxSuccess = 10000}
    (traceHermProp gen :: DensityMatrix -> Property)
  putStrLn "Testing pure state generation:"
  QC.quickCheckWith
    QC.stdArgs {QC.maxSuccess = 10000}
    (traceHermProp gen :: PureStateVector -> Property)

arbPropPure :: PureStateVector -> Property
arbPropPure (PureStateVector sv) = QC.property $ abs (1 - norm_2 sv) < 1e-12

arbPropDensityMatrix :: DensityMatrix -> Property
arbPropDensityMatrix (DensityMatrix dm) =
  QC.property $
  (abs (1 - (magnitude . sumElements . takeDiag) dm) < 1e-12) &&
  (norm_2 (dm - tr dm) < 1e-12)

testStateArb :: IO ()
testStateArb = do
  putStrLn "Testing density matrix arbitrary:"
  QC.quickCheckWith
    QC.stdArgs {QC.maxSuccess = 1000}
    (arbPropDensityMatrix :: DensityMatrix -> Property)
  putStrLn "Testing pure state arbitrary:"
  QC.quickCheckWith
    QC.stdArgs {QC.maxSuccess = 10000}
    (arbPropPure :: PureStateVector -> Property)
