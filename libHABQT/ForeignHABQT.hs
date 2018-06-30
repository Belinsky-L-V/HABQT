module ForeignHABQT where

import Data.Bifunctor (bimap)
import Data.Complex (Complex(..))
import Data.List (unzip)
import qualified Data.Vector as V
import Foreign
import Foreign.C
import HABQTlib.Data
import HABQTlib.MeasurementProcessing (PurePOVM)
import qualified Numeric.LinearAlgebra as LA

tdm :: DensityMatrix
tdm = DensityMatrix $ LA.real $ (2 LA.>< 2) [1 ..]

type CMatrix = [[CDouble]]

convertDM :: DensityMatrix -> (CMatrix, CMatrix)
convertDM (DensityMatrix dm) =
  let (rp, ip) = LA.fromComplex dm
      conv = (map . map) realToFrac . LA.toLists
   in (conv rp, conv ip)

unmarshallSV :: Dim -> Ptr CDouble -> Ptr CDouble -> IO PureStateVector
unmarshallSV dim rPtr iPtr = do
  rl <- peekArray dim rPtr
  il <- peekArray dim iPtr
  let toHM = (dim LA.>< 1) . map realToFrac
  return . PureStateVector $ LA.toComplex (toHM rl, toHM il)

convertWPSV :: WeighedPureStateVector -> ([CDouble], [CDouble])
convertWPSV (WeighedPureStateVector (w, PureStateVector sv)) =
  let sv' = LA.scale (sqrt w :+ 0) sv
      (rp, ip) = LA.fromComplex sv'
      conv = map realToFrac . LA.toList . head . LA.toColumns
   in (conv rp, conv ip)

convertPOVM :: PurePOVM -> ([CDouble], [CDouble])
convertPOVM povm =
  let svs = V.toList povm
      csvs = map convertWPSV svs
   in bimap concat concat . unzip $ csvs

marshallPOVM :: PurePOVM -> Ptr CDouble -> Ptr CDouble -> IO ()
marshallPOVM povm rPtr iPtr = do
  let (r, i) = convertPOVM povm
  pokeArray rPtr r
  pokeArray iPtr i

marshallDM :: DensityMatrix -> Ptr CDouble -> Ptr CDouble -> IO ()
marshallDM dm rPtr iPtr = do
  let (r, i) = convertDM dm
  pokeArray rPtr $ concat r
  pokeArray iPtr $ concat i
