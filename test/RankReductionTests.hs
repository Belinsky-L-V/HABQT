module RankReductionTests
  ( testRankReduction
  ) where

import Data.List (sort)
import HABQTlib.Data
import HABQTlib.Data.Particle
import qualified Numeric.LinearAlgebra as LA
import StateGenTests (arbPropDensityMatrix)
import qualified System.Random.MWC as MWC
import qualified Test.QuickCheck as QC
import Test.QuickCheck
  ( Args(..)
  , Positive(..)
  , Property
  , (.&&.)
  , (===)
  , ioProperty
  , quickCheckWith
  , stdArgs
  )
import TestHelpers

sortSVDprop :: Positive Int -> Positive Int -> Property
sortSVDprop (Positive r) (Positive c) =
  ioProperty $ do
    r1 <- LA.randn r c
    r2 <- LA.randn r c
    let rm = LA.toComplex (r1, r2)
        (u, s) = LA.leftSV rm
        sorted = sort (LA.toList s) == (reverse . LA.toList) s
        normU = (LA.norm_2 . head . LA.toColumns) u
    return $ sorted && abs (normU - 1) < 1e-6

truncateProp' :: Rank -> WeighedDensityMatrix -> Property
truncateProp' rank wdm@(WeighedDensityMatrix (w, _)) =
  let WeighedDensityMatrix (wt, dmt) = truncateRank rank wdm
      preservesWeight = abs (2 * (w - wt) / (w + wt)) <= 1e-12
      (_, s, _) = LA.compactSVD . getDensityMatrix $ dmt
      setsRank = rank === LA.size s
      validDM = arbPropDensityMatrix dmt
   in preservesWeight .&&. setsRank .&&. validDM

truncateProp :: Property
truncateProp =
  let gen = do
        (dim, rank) <- arbDimRank
        wdm <- arbWDM dim rank
        newRank <- QC.suchThat QC.arbitrary (\x -> x > 0 && x <= rank)
        return (newRank, wdm)
   in QC.forAll gen (uncurry truncateProp')

testTruncation :: IO ()
testTruncation = do
  putStrLn "Testing density matrix rank reduction:"
  quickCheckWith stdArgs {maxSuccess = 1000} truncateProp

particleReductionProp :: MWC.GenIO -> Positive Int -> Positive Int -> Property
particleReductionProp gen (Positive dim) (Positive num) =
  ioProperty $ do
    r <- MWC.uniformR (1, dim) gen
    particles <- genParticles dim r num
    let WeighedDensityMatrix (w, dm) = reduceParticlesToMean particles
        dimProp = dim === LA.rows (getDensityMatrix dm)
        rankProp = r === getRank dm
        validProp = arbPropDensityMatrix dm
        weightProp = abs (w - 1) < 1e-12
    return $ dimProp .&&. rankProp .&&. validProp .&&. weightProp

testSVDsort :: IO ()
testSVDsort = do
  putStrLn "Testing singular value sort:"
  quickCheckWith stdArgs {maxSuccess = 1000} sortSVDprop
  putStrLn "Testing particle vector reduction:"
  gen <- MWC.createSystemRandom
  quickCheckWith stdArgs {maxSuccess = 1000} (particleReductionProp gen)

testRankReduction :: IO ()
testRankReduction = do
  testSVDsort
  testTruncation
