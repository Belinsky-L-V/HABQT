module SuperpositionSemigroupTests
  ( testWeighedDensityMatrixSemigroup
  ) where

import HABQTlib.Data
import Test.QuickCheck
  ( Arbitrary(..)
  , Args(..)
  , Gen
  , Property
  , forAll
  , getSize
  , quickCheckWith
  , resize
  , stdArgs
  , suchThat
  )
import TestHelpers

type WDM = WeighedDensityMatrix

testAssocWithGen :: Gen (WDM, WDM, WDM) -> Property
testAssocWithGen gen =
  forAll gen (\(x, y, z) -> ((x <+> y) <+> z) <==> (x <+> (y <+> z)))

gen3DM :: Gen (WDM, WDM, WDM)
gen3DM = do
  n <- suchThat getSize (> 0)
  dm1 <- resize n arbitrary
  dm2 <- resize n arbitrary
  dm3 <- resize n arbitrary
  return (dm1, dm2, dm3)

testWeighedDensityMatrixSemigroup :: IO ()
testWeighedDensityMatrixSemigroup = do
  putStrLn "Testing weighed density matrix semigroup:"
  quickCheckWith stdArgs {maxSuccess = 1000} (testAssocWithGen gen3DM)
