module FidelityTests
  ( testFidelity
  ) where

import HABQTlib.Data
import qualified Test.QuickCheck as QC
import TestHelpers

fidProp :: QC.Positive Dim -> QC.Property
fidProp (QC.Positive dim) =
  let gen = do
        v1 <- QC.resize dim QC.arbitrary
        v2 <- QC.resize dim QC.arbitrary
        return (v1, v2)
      fidMatch :: (PureStateVector, PureStateVector) -> QC.Property
      fidMatch (sv1, sv2) =
        QC.property $
        abs (fidelity sv1 sv2 - fidelityDM (svToDM sv1) (svToDM sv2)) <= 1e-12
   in QC.forAll gen fidMatch

testFidelity :: IO ()
testFidelity = do
  putStrLn "Testing density matrix fidelity:"
  QC.quickCheckWith QC.stdArgs {QC.maxSuccess = 1000} fidProp
