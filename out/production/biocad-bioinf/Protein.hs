{-# LANGUAGE OverloadedStrings #-}

module Protein
    ( putHydrophobicGroups
    ) where

import           Bio.PDB                (modelsFromPDBFile)
import           Bio.Structure          (Model(..), Chain(..), Residue(..), Atom(..))
import           Data.Either            (fromRight)
import           Data.Vector            (Vector)
import           Data.Text              (Text)
import           Linear.V3              (V3(..), _y, _x, _z)
import           Lens.Micro             ((^.))
import           Data.Maybe             (fromMaybe)
import qualified Data.Vector as V
import qualified Data.Text as T

type Radius = Float
type AtomCoords = V3 Float
type Residues = Vector Residue

radius :: Radius
radius = 4

hydrophobicNames :: Vector Text
hydrophobicNames = V.fromList ["LEU", "LIE", "PHE", "TRP", "VAL", "MET", "CYS", "TYR", "ALA"]

findNeighbors :: Radius -> Residues -> Residue -> Residues
findNeighbors r' rs r = V.filter (fromMaybe False . isResiduesNeighbors r' r)  rs

residueCenter :: Residue -> Maybe Atom
residueCenter r = V.find (\a -> atomName a =="CA") (resAtoms r)

showResidue :: Residue -> Text
showResidue res = resName res <> T.pack (show $ resNumber res)

allResiduesFromModels :: Vector Model -> Residues
allResiduesFromModels m = V.concatMap chainResidues $ V.concatMap modelChains m

hydrophobicResiduesFromModels :: Vector Model -> Residues
hydrophobicResiduesFromModels m = V.filter (\v -> V.any (== resName v) hydrophobicNames) (allResiduesFromModels m)

isAtomsNeighbors :: Radius -> Atom -> Atom -> Bool
isAtomsNeighbors r a1 a2 = (x a1 - x a2)**2 + (y a1 - y a2)**2 + (z a1 - z a2)**2 < r**2
  where
    x a = atomCoords a ^. _x
    y a = atomCoords a ^. _y
    z a = atomCoords a ^. _z

isResiduesNeighbors :: Radius -> Residue -> Residue -> Maybe Bool
isResiduesNeighbors r res1 res2 = do
  a1 <- residueCenter res1
  a2 <- residueCenter res2
  return $ isAtomsNeighbors r a1 a2

breadSearch :: Radius -> Residues -> Residues
breadSearch r residues = go V.empty $ V.head residues
  where
    go founded residue =
      if V.elem residue founded
      then founded
      else V.foldl go (V.cons residue founded) (findNeighbors r residues residue)

excludeFounded :: Residues -> Residues -> Residues
excludeFounded x y = V.filter (`V.notElem` y) x

mergeToGroups :: Radius -> Residues -> Vector Residues
mergeToGroups r = go V.empty
  where
    f result founded residues' =
      go (V.cons result founded) $ excludeFounded residues' result
    go founded residues =
      if V.length residues == 0
      then founded
      else f (breadSearch r residues) founded residues

groupName :: Residues -> Text
groupName = V.foldr (\res acc -> acc <> showResidue res) ""

groupCoords :: Residues -> Vector AtomCoords
groupCoords = V.mapMaybe (fmap atomCoords . residueCenter)

putHydrophobicGroups :: IO ()
putHydrophobicGroups = do
  eitherPDB <- modelsFromPDBFile "pbd/1mqk.pdb"
--  eitherPDB <- modelsFromPDBFile "pbd/8cvz.pdb"
--  eitherPDB <- modelsFromPDBFile "pbd/7xk5.pdb"
  let (_, models) = fromRight undefined eitherPDB
      residues = hydrophobicResiduesFromModels models
      resultResidues = mergeToGroups radius residues
      resultNames = V.map groupName resultResidues
--      resultCoords = V.map groupCoords resultResidues
  V.mapM_ print resultNames

