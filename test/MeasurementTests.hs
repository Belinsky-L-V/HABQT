module MeasurementTests
  ( testMeasurements
  ) where

import Data.Maybe (fromJust)
import qualified Data.Vector as V
import HABQTlib.Data
import HABQTlib.MeasurementProcessing
import Numeric.LinearAlgebra (Complex(..), (><))
import qualified Numeric.LinearAlgebra as LA
import StateGenTests (arbPropPure)
import qualified Test.QuickCheck as QC
import Test.QuickCheck (Property, (.&&.))
import TestHelpers

povmValidProp :: PurePOVM -> Property
povmValidProp wsts = normedContents .&&. sumsToUnity
  where
    dim :: Int
    dim =
      (\(WeighedPureStateVector (_, PureStateVector sv)) -> LA.rows sv) $
      V.head wsts
    normedContents =
      QC.conjoin . V.toList $
      V.map (\(WeighedPureStateVector (_, sv)) -> arbPropPure sv) wsts
    msum :: LA.Matrix (LA.Complex Double)
    msum =
      V.foldl'
        (\acc (WeighedPureStateVector (w, PureStateVector sv)) ->
           acc + LA.scalar (w :+ 0) * (sv LA.<> LA.tr sv))
        ((dim >< dim) (repeat (0 :+ 0)))
        wsts
    sumsToUnity = LA.norm_Frob (msum - LA.ident dim) <= 1e-6

povmProbabilityNorm :: PurePOVM -> DensityMatrix -> Property
povmProbabilityNorm povm dm =
  QC.property $ abs (1 - V.sum (measurementProbs povm dm)) <= 1e-8

blochHelper :: Double -> Double -> SingleQbParam
blochHelper = (,)

coordProp :: Property
coordProp =
  let prop sv = fidelity sv sv' > 0.99999
        where
          sv' = blochAnglesToSV . fromJust . svToAngles $ sv
   in QC.forAll (QC.resize 2 QC.arbitrary) prop

testBlochCoords :: IO ()
testBlochCoords = do
  putStrLn "Testing Bloch coordinate transofrmations:"
  QC.quickCheckWith QC.stdArgs {QC.maxSuccess = 1000} coordProp

testAntipodal1QbPOVM :: IO ()
testAntipodal1QbPOVM = do
  putStrLn "Testing antipodal single qubit POVM:"
  QC.quickCheckWith QC.stdArgs {QC.maxSuccess = 1000} $
    ((povmValidProp . mkAntipodalPOVM) .) . blochHelper

genProdPOVM :: Int -> QC.Gen PurePOVM
genProdPOVM qbnum = do
  paramVs' <- QC.vectorOf qbnum QC.arbitrary
  let paramVs = fmap (uncurry blochHelper) paramVs'
      svs = map mkAntipodalPOVM paramVs
  return $ productPOVM svs

testProductPOVM :: IO ()
testProductPOVM = do
  putStrLn "Testing productPOVMs"
  let gen = do
        qbnum <- QC.suchThat QC.arbitrary (\n -> n >= 1 && n <= 6)
        let dim :: Dim
            dim = 2 ^ qbnum
        rank <- QC.suchThat QC.arbitrary (\n -> n <= dim && n > 0)
        povm <- genProdPOVM qbnum
        dm <- arbDM dim rank
        return (povm, dm)
  QC.quickCheckWith QC.stdArgs {QC.maxSuccess = 10000} $
    QC.forAll
      gen
      (\(povm, dm) -> povmValidProp povm .&&. povmProbabilityNorm povm dm)

testMeasurements :: IO ()
testMeasurements = do
  testBlochCoords
  testAntipodal1QbPOVM
  testProductPOVM
