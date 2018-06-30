{-# LANGUAGE ForeignFunctionInterface #-}

module LibHABQT where

import Control.Monad.State.Lazy
import Data.IORef
import Data.Validation
import Foreign
import Foreign.C
import ForeignHABQT
import HABQTlib.Data
import HABQTlib.Data.Particle
import HABQTlib.UnsafeAPI
import qualified System.IO as IO
import qualified System.Random.MWC as MWC

type TomForegin
   = Ptr CDouble -- measurement SV real
      -> Ptr CDouble -- measurement SV imag
          -> Ptr CDouble -- estimate dm real
              -> Ptr CDouble -- estimate dm imag
                  -> Ptr CDouble -- next POVM real
                      -> Ptr CDouble -- next POVM imag
                          -> IO ()

type ArgStorage = (QBitNum, OutputVerb, MHMCiter, OptIter, MWC.GenIO)

type Storage = IORef (ArgStorage, (ParticleHierarchy, [PureStateVector]))

type InitFun = CInt -> CInt -> CInt -> CInt -> CInt -> IO (StablePtr Storage)

foreign export ccall "tomography_init" tomInit :: InitFun

tomInit :: InitFun
tomInit qn' pc' mi' oi' verb' = do
  let inp = validateInputs qn' pc' mi' oi' verb'
  gen <- MWC.createSystemRandom
  case inp of
    Failure errs -> do
      IO.hSetBuffering IO.stderr IO.LineBuffering
      IO.hPutStrLn IO.stderr (unlines errs)
      IO.hPutStrLn
        IO.stderr
        "Free the resulting pointer using tomography_free and re-initialise with valid inputs"
      IO.hFlush IO.stderr
      ph <- initialiseParticleHierarchy 2 1
      mem <- newIORef ((1, NoOutput, 1, 1, gen), (ph, []))
      newStablePtr mem
    Success (qn, pc, mi, oi, verb) -> do
      let dim = 2 ^ qn
          out =
            case verb of
              0 -> NoOutput
              1 -> FullOutput
      ph <- initialiseParticleHierarchy dim pc
      mem <- newIORef ((qn, out, mi, oi, gen), (ph, []))
      newStablePtr mem

foreign export ccall "tomography" foreignTomFun
  :: StablePtr Storage -> TomForegin

foreignTomFun :: StablePtr Storage -> TomForegin
foreignTomFun strPtr svrPtr sviPtr dmrPtr dmiPtr povmrPtr povmiPtr = do
  mem <- deRefStablePtr strPtr
  ((qn, out, mi, oi, gen), s) <- readIORef mem
  let dim = 2 ^ qn
  sv <- unmarshallSV dim svrPtr sviPtr
  case validSV sv of
    Success _ -> do
      ((dm, nextPOVM), sn) <- runStateT (tomographyFun' qn mi oi out gen sv) s
      writeIORef mem ((qn, out, mi, oi, gen), sn)
      marshallDM dm dmrPtr dmiPtr
      marshallPOVM nextPOVM povmrPtr povmiPtr
    Failure [msg] -> do
      IO.hSetBuffering IO.stderr IO.LineBuffering
      IO.hPutStrLn IO.stderr $ msg ++ " Doing nothing."
      IO.hFlush IO.stderr

foreign export ccall "tomography_free" tomFree :: StablePtr Storage -> IO ()

tomFree :: StablePtr Storage -> IO ()
tomFree = freeStablePtr

validateInputs ::
     CInt
  -> CInt
  -> CInt
  -> CInt
  -> CInt
  -> Validation [String] (QBitNum, NumberOfParticles, MHMCiter, OptIter, Int)
validateInputs qn' pn' mi' oi' verb' =
  let qn = fromIntegral qn'
      pn = fromIntegral pn'
      verb = fromIntegral verb'
      mi = fromIntegral mi'
      oi = fromIntegral oi'
      qnM = ["Number of quantum bits must be a positive integer."]
      pnM = ["Number of particles per rank must be a positive integer."]
      vM =
        [ "Verbosity level must be an integer equal to 0 (no output) or 1 (full output)."
        ]
      miM = ["Number of MHMC iterations must be a positive integer."]
      vQn = validate qnM (> 0) qn
      vPn = validate pnM (> 0) pn
      vVerb = validate vM ((||) <$> (== 0) <*> (== 1)) verb
      vMi = validate miM (> 0) mi
      vOi = validOptIter oi
   in (,,,,) <$> vQn <*> vPn <*> vMi <*> vOi <*> vVerb
