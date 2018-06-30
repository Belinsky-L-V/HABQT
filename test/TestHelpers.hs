module TestHelpers where

import Control.Newtype.Generics (over)
import Data.Complex
import qualified Data.Vector as V
import HABQTlib.Data
import HABQTlib.Data.Particle
import qualified Numeric.LinearAlgebra as LA
import qualified Test.QuickCheck as QC
import Test.QuickCheck (Arbitrary(..), Gen, Positive(..))

getRank :: DensityMatrix -> Rank
getRank (DensityMatrix dm) =
  let (_, s, _) = LA.compactSVD dm
   in LA.size s

arbDimRank :: QC.Gen (Dim, Rank)
arbDimRank = do
  dim <- QC.suchThat arbitrary (> 1)
  Positive rank <- QC.suchThat arbitrary (\(Positive v) -> v <= dim)
  return (dim, rank)

arbDM :: Dim -> Rank -> QC.Gen DensityMatrix
arbDM dim rank = do
  let genXs =
        QC.vectorOf
          dim
          (QC.vectorOf rank (QC.arbitrary :: QC.Gen (Complex Double)))
  xs <- QC.suchThat genXs (\ll -> any (/= 0.0 :+ 0.0) (zipWith (!!) ll [0 ..]))
  let a = LA.fromLists xs
      h = a LA.<> LA.tr a
      hTr = LA.sumElements $ LA.takeDiag h
      dm = h / LA.scalar hTr
  return $ DensityMatrix dm

arbWDM :: Dim -> Rank -> QC.Gen WeighedDensityMatrix
arbWDM dim rank = do
  Positive w <- QC.arbitrary
  dm <- arbDM dim rank
  return $ WeighedDensityMatrix (w, dm)

arbNormedVecDim :: Int -> Gen (LA.Matrix (Complex Double))
arbNormedVecDim dim = do
  let genXs = QC.vectorOf dim (arbitrary :: Gen (Complex Double))
  xs <- QC.suchThat genXs (any (/= 0))
  let sv = LA.asColumn . LA.fromList $ xs
  return $ sv / LA.toComplex (LA.scalar (LA.norm_2 sv), 0)

arbitraryWeighed :: Arbitrary x => Gen (Weight, x)
arbitraryWeighed = do
  dim <- QC.suchThat QC.getSize (> 0)
  w <- getPositive <$> arbitrary
  dm <- QC.resize dim arbitrary
  return (w, dm)

arbParticles :: Dim -> Rank -> NumberOfParticles -> QC.Gen Particles
arbParticles dim rank num = do
  vdms <- V.replicateM num (arbWDM dim rank)
  let wr = V.foldl' (\acc (WeighedDensityMatrix (wi, _)) -> acc + wi) 0 vdms
      vdmsn = V.map (over WeighedDensityMatrix (\(w, dm) -> (w / wr, dm))) vdms
  QC.Positive w <- QC.arbitrary
  return $ Particles rank w num vdmsn

instance Arbitrary DensityMatrix where
  arbitrary = do
    dim <- QC.suchThat QC.getSize (> 0)
    arbDM dim dim

instance Arbitrary PureStateVector where
  arbitrary = do
    dim <- QC.suchThat QC.getSize (> 0)
    sv <- arbNormedVecDim dim
    return . PureStateVector $ sv

instance Arbitrary WeighedDensityMatrix where
  arbitrary = do
    (w, dm) <- arbitraryWeighed
    return . WeighedDensityMatrix $ (w, dm)

instance Arbitrary WeighedPureStateVector where
  arbitrary = do
    (w, sv) <- arbitraryWeighed
    return . WeighedPureStateVector $ (w, sv)

instance QC.Arbitrary Particles where
  arbitrary = do
    (dim, rank) <- arbDimRank
    num <- QC.suchThat QC.arbitrary (> 100)
    arbParticles dim rank num

class MEq a where
  infix 4 <==>
  (<==>) :: a -> a -> Bool

instance MEq DensityMatrix where
  DensityMatrix dm1 <==> DensityMatrix dm2 = LA.norm_Frob (dm1 - dm2) < 1e-6

instance MEq WeighedDensityMatrix where
  WeighedDensityMatrix (w1, DensityMatrix dm1) <==> WeighedDensityMatrix (w2, DensityMatrix dm2) =
    LA.norm_Frob (LA.scalar (w1 :+ 0) * dm1 - LA.scalar (w2 :+ 0) * dm2) < 1e-6

instance MEq PureStateVector where
  sv1 <==> sv2 = 1 - fidelity sv1 sv2 < 1e-6

instance MEq WeighedPureStateVector where
  WeighedPureStateVector (w1, sv1) <==> WeighedPureStateVector (w2, sv2) =
    abs (w1 - w2) <= 1e-6 * 0.5 * (w1 + w2) && sv1 <==> sv2
