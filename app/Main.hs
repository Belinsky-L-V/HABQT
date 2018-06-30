{-# LANGUAGE RecordWildCards #-}

module Main where

import Data.Semigroup ((<>))
import HABQTlib.Data
import HABQTlib.UnsafeAPI
import Options.Applicative
import Streaming (Of, Stream)
import qualified Streaming.Prelude as S
import qualified System.IO as IO

data CLIargs = CLIargs
  { cliQbNum :: QBitNum
  , cliRank :: Rank
  , cliMeasNum :: Int
  , cliExpNum :: Int
  , cliPtNum :: NumberOfParticles
  , cliVerb :: Int
  , cliFilePath :: IO.FilePath
  , cliMHMCiter :: MHMCiter
  , cliOptIter :: OptIter
  }

msgPositive :: String
msgPositive = "must be a positive integer"

readInt :: (Int -> Bool) -> String -> String -> Either String Int
readInt cond msg s =
  let pstr :: [(Int, String)]
      pstr = reads s
      go [(i, "")] =
        if cond i
          then Right i
          else Left msg
      go _ = Left msg
   in go pstr

readPostiveInt :: String -> Either String Int
readPostiveInt = readInt (> 0) msgPositive

positive :: ReadM Int
positive = eitherReader readPostiveInt

cliparse :: Parser CLIargs
cliparse =
  CLIargs <$>
  option
    positive
    (long "quantum-bit-number" <> short 'q' <>
     help "Number of quantum bits under tomography" <>
     showDefault <>
     value 1 <>
     metavar "QBNUM") <*>
  option
    positive
    (long "rank" <> short 'r' <> help "Rank of true states" <> showDefault <>
     value 1 <>
     metavar "RANK") <*>
  option
    positive
    (long "measurements-per-experiment" <> short 'm' <>
     help "Number of measurements to make per-experiment" <>
     showDefault <>
     value 100 <>
     metavar "MNUM") <*>
  option
    positive
    (long "experiment-number" <> short 'e' <>
     help "Number of tomography experiments to run" <>
     showDefault <>
     value 10 <>
     metavar "EXPNUM") <*>
  option
    positive
    (long "particle-number" <> short 'p' <>
     help
       "Number of particles to use for approximating the distribution over states (per rank)" <>
     showDefault <>
     value 1000 <>
     metavar "PNUM") <*>
  option
    auto
    (long "verbosity" <> short 'v' <>
     help
       "Verbosity level of output to stdout from 0 (no output) to 2 (full output)" <>
     showDefault <>
     value 2 <>
     metavar "VERB") <*>
  strOption
    (long "output-file-path" <> short 'o' <>
     help
       "Path to file which infidelities of tomographic estimates will be appended to" <>
     showDefault <>
     value "./output.txt" <>
     metavar "OUTPATH") <*>
  option
    positive
    (long "mhmc-iterations" <>
     help
       "Number of Metropolis–Hastings steps to perform when resampling (after adjusting for acceptance rate)" <>
     showDefault <>
     value 50 <>
     metavar "MHMCITER") <*>
  option
    positive
    (long "opt-iterations" <>
     help
       "Number of optimisation steps to perform when searching for optimal measurment" <>
     showDefault <>
     value 50 <>
     metavar "POVMITER")

writeCsInfid ::
     Int
  -> Int
  -> QBitNum
  -> Rank
  -> NumberOfParticles
  -> MHMCiter
  -> OptIter
  -> OutputVerb
  -> IO.FilePath
  -> IO ()
writeCsInfid numRuns numMeas qbn rank pn mi oi v fp = do
  let tomStr = streamTo numMeas (streamResults' qbn rank pn mi oi v)
      multStr = S.replicateM numRuns . tomStr
      dupIO fh s = do
        IO.hPutStr fh s
        IO.hFlush IO.stdout
  IO.withFile fp IO.AppendMode (S.effects . S.effects . multStr . dupIO)

streamTo ::
     Int
  -> Stream (Of Double) IO ()
  -> (String -> IO ())
  -> Stream (Of ()) IO ()
streamTo n is sf = do
  let rs = S.map show is
      fs = S.take 1 rs
      rest = S.drop 1 rs
  S.mapM sf . S.yield $ "\n"
  S.mapM sf fs
  S.mapM (\x -> sf (", " ++ x)) . S.take n $ rest
  S.mapM sf . S.yield $ "\n"

main :: IO ()
main = do
  CLIargs {..} <- execParser opts
  r' <-
    if cliRank > 2 ^ cliQbNum
      then do
        putStrLn "WARNING: RANK > 2 ^ QBNUM, truncating to 2 ^ QBNUM"
        return (2 ^ cliQbNum)
      else return cliRank
  v <-
    case cliVerb of
      0 -> return NoOutput
      1 -> return FidOutput
      2 -> return FullOutput
      _ -> do
        putStrLn
          "WARNING: VERB isn't equal to 0,1 or 2, defaulting to full output"
        return FullOutput
  writeCsInfid
    cliExpNum
    cliMeasNum
    cliQbNum
    r'
    cliPtNum
    cliMHMCiter
    cliOptIter
    v
    cliFilePath
  where
    opts =
      info
        (cliparse <**> helper)
        (fullDesc <>
         progDesc "Simulate HABQT and write infidelities of estimates to file" <>
         header "Hierarchical Adaptive Bayesian Quantum Tomography simulation")
