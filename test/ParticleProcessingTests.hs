{-# LANGUAGE NamedFieldPuns #-}

module ParticleProcessingTests
  ( testParticleHandling
  ) where

import Control.Applicative (liftA2)
import Data.Complex (Complex(..))
import Data.List (sort)
import qualified Data.Vector as V
import HABQTlib.Data
import HABQTlib.Data.Particle
import HABQTlib.RandomStates
import qualified Numeric.LinearAlgebra as LA
import StateGenTests (arbPropDensityMatrix)
import qualified Streaming.Prelude as S
import qualified System.Random.MWC as MWC
import qualified Test.QuickCheck as QC
import Test.QuickCheck
  ( Args(..)
  , Positive(..)
  , Property
  , (.&&.)
  , arbitrary
  , quickCheckWith
  , stdArgs
  , suchThat
  )
import TestHelpers

lik :: Double -> Double
lik p = p / (1 - p)

loglik :: Double -> Double
loglik = log . lik

likProp :: PureStateVector -> DensityMatrix -> Property
likProp sv dm =
  let l = pureStateLikelihood sv dm
   in QC.property (l >= 0) .&&. QC.property (l <= 1)

testLikelihood :: IO ()
testLikelihood = do
  putStrLn "Testing likelihood computation for density matrices:"
  quickCheckWith stdArgs {maxSuccess = 1000} likProp

expProp :: QC.Blind Particles -> Property
expProp (QC.Blind pts) =
  let dm = snd . getWDM . V.foldl1' (<+>) . ptsParticles $ pts
      dim = LA.rows . getDensityMatrix $ dm
      z = (dim LA.>< dim) $ repeat 0
      scaleDM :: Weight -> DensityMatrix -> LA.Matrix (Complex Double)
      scaleDM w (DensityMatrix dm'') = LA.scale (w :+ 0) dm''
      dm' = DensityMatrix $ foldOverPts id scaleDM (+) z pts
      matchProp = QC.property $ dm <==> dm'
      validProp = arbPropDensityMatrix dm'
   in matchProp .&&. validProp

testExpectation :: IO ()
testExpectation = do
  putStrLn "Testing particle expectation function:"
  quickCheckWith stdArgs {maxSuccess = 500} expProp

particleUpdateFidR :: PureStateVector -> Particles -> Double
particleUpdateFidR testState@(PureStateVector sv) particles =
  let testPure = DensityMatrix $ sv LA.<> LA.tr sv
      originalEstimate = snd . getWDM $ reduceParticlesToMean particles
      originalFid = fidelityDM testPure originalEstimate
      updatedParticles = updateParticles testState particles
      updatedEstimate = snd . getWDM $ reduceParticlesToMean updatedParticles
      updatedFid = fidelityDM testPure updatedEstimate
   in loglik updatedFid - loglik originalFid

updateProp :: Property
updateProp =
  let gen = do
        (dim, rank) <- arbDimRank
        Positive num <- arbitrary :: QC.Gen (Positive NumberOfParticles)
        stateVector <- QC.resize dim arbitrary
        particles <- arbParticles dim rank (fromIntegral num)
        return (stateVector, particles)
      prop (sv, pts) =
        let wtsf Particles {ptsParticles} =
              V.foldl1' (+) . V.map (fst . getWDM) $ ptsParticles
            ws = wtsf pts
            pts' = updateParticles sv pts
            ws' = wtsf pts'
            wsProp = abs (ws - 1) < 1e-10
            ws'Prop = abs (ws' - 1) < 1e-10
         in wsProp .&&. ws'Prop
      msg = "of updates move the estimate closer"
   in QC.forAll
        gen
        (\(sv, pts) ->
           QC.classify (particleUpdateFidR sv pts >= 0) msg (prop (sv, pts)))

ecdfProp :: Particles -> Property
ecdfProp Particles {ptsParticles} =
  let cdf = ecdf ptsParticles
      hp = V.head cdf > 0
      tp = abs (V.last cdf - 1) < 1e-10
      cdfSorted = V.fromList . sort . V.toList $ cdf
      sorted = cdf == cdfSorted
   in hp .&&. tp .&&. sorted

icdfProp :: Particles -> Property
icdfProp Particles {ptsParticles} =
  QC.ioProperty $ do
    gen <- MWC.createSystemRandom
    x <- MWC.uniform gen
    let cdf = ecdf ptsParticles
        idx = icdf cdf x
        Just idx' = V.findIndex (> x) cdf
    return (idx == idx')

hDistProp :: ParticleHierarchy -> PureStateVector -> Double
hDistProp ph sv =
  let dm = DensityMatrix $ getStateVector sv LA.<> LA.tr (getStateVector sv)
      originalEstimate = getMixedEstimate ph
      updatedEstimate = getMixedEstimate $ updateParticleHierarchy sv ph
      originalFid = fidelityDM dm originalEstimate
      updatedFid = fidelityDM dm updatedEstimate
   in loglik updatedFid - loglik originalFid

hierarchyBatchProp ::
     Int -> Positive Dim -> Positive NumberOfParticles -> IO Property
hierarchyBatchProp batchSize (Positive dim) (Positive num) = do
  ph <- initialiseParticleHierarchy dim num
  passed <-
    S.length_ . S.filter (>= 0) . S.replicateM batchSize $
    hDistProp ph <$> genPureSV dim
  return . QC.property $
    fromIntegral passed / (fromIntegral batchSize :: Double) >= 0.95

hierarchyUpdateProp :: Int -> Property
hierarchyUpdateProp batchSize =
  let gen = do
        dim <- suchThat arbitrary (\x -> x > 1 && x < 20)
        num <- suchThat arbitrary (> 100)
        return (Positive dim, Positive num)
   in QC.forAll
        gen
        (\(dim, num) -> QC.ioProperty $ hierarchyBatchProp batchSize dim num)

effectiveSampleSizeProp :: Property
effectiveSampleSizeProp =
  let genP = do
        (dim, rank) <- arbDimRank
        Positive num <- arbitrary :: QC.Gen (Positive NumberOfParticles)
        arbParticles dim rank (fromIntegral num)
      sizeP ps@(Particles _ _ n _) =
        QC.classify
          (es <= fromIntegral n / 2)
          "of effective sizes are less than half of total particle number"
          ((es <= fromIntegral n) && (es >= 0))
        where
          es = effectiveSize ps
   in QC.forAll genP sizeP

nudgeProp :: Property
nudgeProp =
  let gen = do
        (dim, rank) <- arbDimRank
        dm <- arbDM dim rank
        return (dim, rank, dm)
      nudged (dim, rank, dm) =
        truncateRank rank <$> nudgeParticle dim 1e-2 (mkWDM1 dm)
      prop (dim, rank, dm) =
        QC.ioProperty $ do
          WeighedDensityMatrix (w, ndm) <- nudged (dim, rank, dm)
          let validDM = arbPropDensityMatrix ndm
              preservesRank = getRank dm == getRank ndm
              preservesWeight = w == 1
          return $
            QC.classify
              (fidelityDM ndm dm > 0.99)
              "fidelities exceed 99%"
              (validDM .&&. preservesRank .&&. preservesWeight)
   in QC.forAll gen prop

type ResamplingFunIO = MWC.GenIO -> Particles -> IO Particles

resampleProp :: MWC.GenIO -> ResamplingFunIO -> Property
resampleProp gen resamplingFun =
  let genA = do
        dim <- suchThat (arbitrary :: QC.Gen Dim) (liftA2 (&&) (< 10) (> 1))
        Positive rank <-
          suchThat (arbitrary :: QC.Gen (Positive Rank)) (< Positive dim)
        num <- suchThat arbitrary (> 1000)
        QC.Blind <$> arbParticles dim rank num
      preservesMeanIO :: QC.Blind Particles -> IO Property
      preservesMeanIO (QC.Blind pts) = do
        rpts <- resamplingFun gen pts
        let WeighedDensityMatrix (w1, dm1) = reduceParticlesToMean pts
            WeighedDensityMatrix (w2, dm2) = reduceParticlesToMean rpts
            sumw wacc (WeighedDensityMatrix (w', _)) = wacc + w'
            getpw Particles {ptsParticles = ps} = V.foldl' sumw 0 ps
            w1' = getpw pts
            w2' = getpw rpts
            propWeight = QC.property $ 2 * abs (w1 - w2) / (w1 + w2) < 1e-12
            propWeight' =
              QC.property $ 2 * abs (w1' - w2') / (w1' + w2') < 1e-12
            Particles _ _ num _ = rpts
            propESabs =
              QC.property $ abs (effectiveSize rpts - fromIntegral num) < 1
            propESincrease =
              QC.property $ effectiveSize rpts > effectiveSize pts
            propAll =
              propWeight .&&. propWeight' .&&. propESincrease .&&. propESabs
        return $
          QC.classify
            (fidelityDM dm1 dm2 > 0.99)
            "fidelities between pre- and post-resampling means exceed 99%"
            propAll
   in QC.forAll genA (QC.ioProperty . preservesMeanIO)

testParticleResampling :: IO ()
testParticleResampling = do
  gen <- MWC.createSystemRandom
  putStrLn "Testing multinomial particle resampling without checks:"
  quickCheckWith stdArgs {maxSuccess = 100} $ resampleProp gen resampleMultinom

testParticleNudging :: IO ()
testParticleNudging = do
  putStrLn "Testing particle nudging:"
  quickCheckWith stdArgs {maxSuccess = 1000} nudgeProp

testEffectiveSampleSize :: IO ()
testEffectiveSampleSize = do
  putStrLn "Testing calculation of effective sample size:"
  quickCheckWith stdArgs {maxSuccess = 1000} effectiveSampleSizeProp

testParticleUpdate :: IO ()
testParticleUpdate = do
  putStrLn "Testing update of particles:"
  quickCheckWith stdArgs {maxSuccess = 1000} updateProp

testEcdf :: IO ()
testEcdf = do
  putStrLn "Testing ecdf calculation:"
  quickCheckWith stdArgs {maxSuccess = 10000} ecdfProp
  putStrLn "Testing iecdf calculation:"
  quickCheckWith stdArgs {maxSuccess = 10000} icdfProp

testHierarchyUpdate :: IO ()
testHierarchyUpdate = do
  putStrLn "Testing update of particle hierarchies:"
  quickCheckWith stdArgs {maxSuccess = 50} (hierarchyUpdateProp 100)

testParticleHandling :: IO ()
testParticleHandling = do
  testLikelihood
  testExpectation
  testParticleUpdate
  testEcdf
  testHierarchyUpdate
  testEffectiveSampleSize
  testParticleNudging
  testParticleResampling
