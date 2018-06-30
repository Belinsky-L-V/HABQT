import FidelityTests
import MeasurementTests
import ParticleProcessingTests
import RankReductionTests
import StateGenTests
import SuperpositionSemigroupTests

main :: IO ()
main = do
  putStrLn ""
  testStateGen
  testStateArb
  testFidelity
  testWeighedDensityMatrixSemigroup
  testRankReduction
  testParticleHandling
  testMeasurements
