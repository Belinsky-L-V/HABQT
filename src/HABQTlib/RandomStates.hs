{-|
Module      : HABQTlib.RandomStates

Generation of random pure (from Haar measure) and mixed states (from measure
induced by partial tracing of purified states).
-}
module HABQTlib.RandomStates
  ( genPureSV
  , genDM
  ) where

import HABQTlib.Data
import Numeric.LinearAlgebra
  ( norm_2
  , randn
  , scalar
  , sumElements
  , takeDiag
  , toComplex
  , tr
  )
import qualified Numeric.LinearAlgebra as LA

-- | Generate a random mixed state of specified rank.
genDM :: Dim -> Rank -> IO DensityMatrix
genDM dim r = do
  r1 <- randn dim r
  r2 <- randn dim r
  let a = toComplex (r1, r2)
  let h = a LA.<> tr a
  let hTr = sumElements $ takeDiag h
  return . DensityMatrix $ h / scalar hTr

-- | Generate a random pure state from Hilbert space of given dimension.
genPureSV :: Dim -> IO PureStateVector
genPureSV dim = do
  r1 <- randn dim 1
  r2 <- randn dim 1
  let sv = toComplex (r1, r2)
  return . PureStateVector $ sv / toComplex (scalar (norm_2 sv), 0)
