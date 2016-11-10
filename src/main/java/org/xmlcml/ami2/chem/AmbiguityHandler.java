package org.xmlcml.ami2.chem;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.xmlcml.ami2.chem.Joinable.JoinPoint;
import org.xmlcml.euclid.Angle;
import org.xmlcml.euclid.Angle.Units;
import org.xmlcml.graphics.svg.SVGLine;

import com.google.common.collect.UnionFind;

public class AmbiguityHandler {
	
	static class MutuallyExclusiveShortLineTriple {
	
		HatchedBond hatchedBond;
		Charge minus;
	
		SVGLine line;
		SingleBond singleBond;
		
		public MutuallyExclusiveShortLineTriple(HatchedBond hatchedBond, Charge minus, SVGLine line) {
			this.hatchedBond = hatchedBond;
			this.minus = minus;
			this.line = line;
		}
	
	}

	static class MutuallyExclusiveShortLinePairTriple {
		
		HatchedBond hatchedBond;
	
		SVGLine line1;
		SVGLine line2;
		DoubleBond doubleBond;
		SingleBond singleBond1;
		SingleBond singleBond2;
		
		public MutuallyExclusiveShortLinePairTriple(HatchedBond hatchedBond, SVGLine line1, SVGLine line2) {
			this.hatchedBond = hatchedBond;
			this.line1 = line1;
			this.line2 = line2;
		}
		
	}

	static class MutuallyExclusiveShortLineTripleTriple {
		
		HatchedBond hatchedBond;
	
		SVGLine line1;
		SVGLine line2;
		SVGLine line3;
		TripleBond tripleBond;
		SingleBond singleBond1;
		SingleBond singleBond2;
		SingleBond singleBond3;
		
		public MutuallyExclusiveShortLineTripleTriple(HatchedBond hatchedBond, SVGLine line1, SVGLine line2, SVGLine line3) {
			this.hatchedBond = hatchedBond;
			this.line1 = line1;
			this.line2 = line2;
			this.line3 = line3;
		}
		
		boolean setIfEquals(TripleBond tripleBond) {
			HashSet<SVGLine> set = new HashSet<SVGLine>(Arrays.asList(line1, line2, line3));
			set.removeAll(tripleBond.getLines());
			if (set.size() == 0) {
				this.tripleBond = tripleBond;
				return true;
			}
			return false;
		}
		
	}
	
	static abstract class LineGroupPair {
		
		abstract UnsaturatedBond getUnsaturatedBond();
		
		abstract List<SingleBond> getSingleBonds();

		abstract List<SVGLine> getLines();
		
	}

	static class MutuallyExclusiveLinePairPair extends LineGroupPair {
	
		SVGLine line1;
		SVGLine line2;
		DoubleBond doubleBond;
		SingleBond singleBond1;
		SingleBond singleBond2;
		
		public MutuallyExclusiveLinePairPair(DoubleBond doubleBond) {
			this.line1 = doubleBond.getLine(0);
			this.line2 = doubleBond.getLine(1);
			this.doubleBond = doubleBond;
		}

		@Override
		UnsaturatedBond getUnsaturatedBond() {
			return doubleBond;
		}

		@Override
		List<SingleBond> getSingleBonds() {
			return Arrays.asList(singleBond1, singleBond2);
		}

		@Override
		List<SVGLine> getLines() {
			return Arrays.asList(line1, line2);
		}
		
	}

	static class MutuallyExclusiveLineTriplePair extends LineGroupPair {
	
		SVGLine line1;
		SVGLine line2;
		SVGLine line3;
		TripleBond tripleBond;
		SingleBond singleBond1;
		SingleBond singleBond2;
		SingleBond singleBond3;
		
		public MutuallyExclusiveLineTriplePair(TripleBond tripleBond) {
			this.line1 = tripleBond.getLine(0);
			this.line2 = tripleBond.getLine(1);
			this.line3 = tripleBond.getLine(2);
			this.tripleBond = tripleBond;
		}

		@Override
		UnsaturatedBond getUnsaturatedBond() {
			return tripleBond;
		}

		@Override
		List<SingleBond> getSingleBonds() {
			return Arrays.asList(singleBond1, singleBond2, singleBond3);
		}

		@Override
		List<SVGLine> getLines() {
			return Arrays.asList(line1, line2, line3);
		}
		
	}

	public List<MutuallyExclusiveShortLineTriple> mutuallyExclusiveShortLineTriples;
	public List<MutuallyExclusiveShortLinePairTriple> mutuallyExclusiveShortLinePairTriples;
	public List<MutuallyExclusiveShortLineTripleTriple> mutuallyExclusiveShortLineTripleTriples;
	public List<MutuallyExclusiveLinePairPair> mutuallyExclusiveLinePairPairs;
	public List<MutuallyExclusiveLineTriplePair> mutuallyExclusiveLineTriplePairs;
	
	private ChemistryBuilder chemistryBuilder;

	public AmbiguityHandler(ChemistryBuilder chemistryBuilder) {
		this.chemistryBuilder = chemistryBuilder;
	}

	private void removeJoinable(Joinable joinable) {
		chemistryBuilder.getHigherPrimitives().getJoinableList().remove(joinable);
		chemistryBuilder.getHigherPrimitives().getDoubleBondList().remove(joinable);
		chemistryBuilder.getHigherPrimitives().getHatchedBondList().remove(joinable);
		chemistryBuilder.getHigherPrimitives().getLineChargeList().remove(joinable);
	}

	void handleAmbiguities(UnionFind<JoinPoint> joinPointsGroupedIntoJunctions) {
		for (MutuallyExclusiveShortLineTriple triple : mutuallyExclusiveShortLineTriples) {
			handleMutuallyExclusiveShortLineTriple(joinPointsGroupedIntoJunctions, triple);
		}
		
		for (MutuallyExclusiveShortLinePairTriple triple : mutuallyExclusiveShortLinePairTriples) {
			handleMutuallyExclusiveShortLinePairTriple(joinPointsGroupedIntoJunctions, triple);
		}
		
		for (MutuallyExclusiveShortLineTripleTriple triple : mutuallyExclusiveShortLineTripleTriples) {
			handleMutuallyExclusiveShortLineTripleTriple(joinPointsGroupedIntoJunctions, triple);
		}

		Set<SingleBond> singleBonds = new HashSet<SingleBond>();
		for (MutuallyExclusiveLinePairPair pair : mutuallyExclusiveLinePairPairs) {
			singleBonds.addAll(pair.getSingleBonds());
		}
		for (MutuallyExclusiveLineTriplePair pair : mutuallyExclusiveLineTriplePairs) {
			singleBonds.addAll(pair.getSingleBonds());
		}
		
		for (MutuallyExclusiveLinePairPair pair : mutuallyExclusiveLinePairPairs) {
			handleMutuallyExclusiveLineGroupPair(joinPointsGroupedIntoJunctions, singleBonds, pair);
		}
		
		for (MutuallyExclusiveLineTriplePair pair : mutuallyExclusiveLineTriplePairs) {
			handleMutuallyExclusiveLineGroupPair(joinPointsGroupedIntoJunctions, singleBonds, pair);
		}
	}

	private void handleMutuallyExclusiveShortLineTriple(UnionFind<JoinPoint> joinPointsGroupedIntoJunctions, MutuallyExclusiveShortLineTriple triple) {
		JoinPoint firstSingleBond = triple.singleBond.getJoinPoints().get(0);
		JoinPoint secondSingleBond = triple.singleBond.getJoinPoints().get(1);
		JoinPoint firstHatchedBond = triple.hatchedBond.getJoinPoints().get(0);
		JoinPoint secondHatchedBond = triple.hatchedBond.getJoinPoints().get(1);
		JoinPoint minus = (triple.minus == null ? null : triple.minus.getJoinPoints().get(0));
		if (joinPointsGroupedIntoJunctions.getSizeOfPartition(firstSingleBond) == 1 && joinPointsGroupedIntoJunctions.getSizeOfPartition(secondSingleBond) == 1 && joinPointsGroupedIntoJunctions.getSizeOfPartition(firstHatchedBond) > 1 && joinPointsGroupedIntoJunctions.getSizeOfPartition(secondHatchedBond) > 1) {
			undoDamageFromIncorrectMinus(joinPointsGroupedIntoJunctions, minus);
			joinPointsGroupedIntoJunctions.remove(firstSingleBond);
			joinPointsGroupedIntoJunctions.remove(secondSingleBond);
			joinPointsGroupedIntoJunctions.remove(minus);
			removeJoinable(triple.singleBond);
			chemistryBuilder.getHigherPrimitives().getLineList().remove(triple.line);
			removeJoinable(triple.minus);
		} else if (joinPointsGroupedIntoJunctions.getSizeOfPartition(firstSingleBond) == 1 && joinPointsGroupedIntoJunctions.getSizeOfPartition(secondSingleBond) == 1 && (joinPointsGroupedIntoJunctions.getSizeOfPartition(firstHatchedBond) == 1 || joinPointsGroupedIntoJunctions.getSizeOfPartition(secondHatchedBond) == 1) && minus != null) {
			joinPointsGroupedIntoJunctions.remove(firstSingleBond);
			joinPointsGroupedIntoJunctions.remove(secondSingleBond);
			joinPointsGroupedIntoJunctions.remove(firstHatchedBond);
			joinPointsGroupedIntoJunctions.remove(secondHatchedBond);
			removeJoinable(triple.singleBond);
			chemistryBuilder.getHigherPrimitives().getLineList().remove(triple.line);
			removeJoinable(triple.hatchedBond);
			removeJoinable(triple.hatchedBond);
		} else {
			undoDamageFromIncorrectMinus(joinPointsGroupedIntoJunctions, minus);
			joinPointsGroupedIntoJunctions.remove(firstHatchedBond);
			joinPointsGroupedIntoJunctions.remove(secondHatchedBond);
			joinPointsGroupedIntoJunctions.remove(minus);
			removeJoinable(triple.hatchedBond);
			removeJoinable(triple.minus);
			removeJoinable(triple.hatchedBond);
			removeJoinable(triple.minus);
		}
	}
	
	private void handleMutuallyExclusiveShortLinePairTriple(UnionFind<JoinPoint> joinPointsGroupedIntoJunctions, MutuallyExclusiveShortLinePairTriple triple) {
		joinPointsGroupedIntoJunctions.remove(triple.singleBond1.getJoinPoints().get(0));
		joinPointsGroupedIntoJunctions.remove(triple.singleBond1.getJoinPoints().get(1));
		joinPointsGroupedIntoJunctions.remove(triple.singleBond2.getJoinPoints().get(0));
		joinPointsGroupedIntoJunctions.remove(triple.singleBond2.getJoinPoints().get(1));
		removeJoinable(triple.singleBond1);
		removeJoinable(triple.singleBond2);
		chemistryBuilder.getHigherPrimitives().getLineList().remove(triple.line1);
		chemistryBuilder.getHigherPrimitives().getLineList().remove(triple.line2);
		if (triple.doubleBond == null) {
			return;
		}
		JoinPoint firstDoubleBond = triple.doubleBond.getJoinPoints().get(0);
		JoinPoint secondDoubleBond = triple.doubleBond.getJoinPoints().get(1);
		JoinPoint firstHatchedBond = triple.hatchedBond.getJoinPoints().get(0);
		JoinPoint secondHatchedBond = triple.hatchedBond.getJoinPoints().get(1);
		if ((joinPointsGroupedIntoJunctions.getSizeOfPartition(firstHatchedBond) > 1 || joinPointsGroupedIntoJunctions.getSizeOfPartition(secondHatchedBond) > 1) && joinPointsGroupedIntoJunctions.getSizeOfPartition(firstDoubleBond) == 1 && joinPointsGroupedIntoJunctions.getSizeOfPartition(secondDoubleBond) == 1) {
			joinPointsGroupedIntoJunctions.remove(firstDoubleBond);
			joinPointsGroupedIntoJunctions.remove(secondDoubleBond);
			removeJoinable(triple.doubleBond);
		} else {
			joinPointsGroupedIntoJunctions.remove(firstHatchedBond);
			joinPointsGroupedIntoJunctions.remove(secondHatchedBond);
			removeJoinable(triple.hatchedBond);
		}
	}
	
	private void handleMutuallyExclusiveShortLineTripleTriple(UnionFind<JoinPoint> joinPointsGroupedIntoJunctions, MutuallyExclusiveShortLineTripleTriple triple) {
		joinPointsGroupedIntoJunctions.remove(triple.singleBond1.getJoinPoints().get(0));
		joinPointsGroupedIntoJunctions.remove(triple.singleBond1.getJoinPoints().get(1));
		joinPointsGroupedIntoJunctions.remove(triple.singleBond2.getJoinPoints().get(0));
		joinPointsGroupedIntoJunctions.remove(triple.singleBond2.getJoinPoints().get(1));
		joinPointsGroupedIntoJunctions.remove(triple.singleBond3.getJoinPoints().get(0));
		joinPointsGroupedIntoJunctions.remove(triple.singleBond3.getJoinPoints().get(1));
		removeJoinable(triple.singleBond1);
		removeJoinable(triple.singleBond2);
		removeJoinable(triple.singleBond3);
		chemistryBuilder.getHigherPrimitives().getLineList().remove(triple.line1);
		chemistryBuilder.getHigherPrimitives().getLineList().remove(triple.line2);
		chemistryBuilder.getHigherPrimitives().getLineList().remove(triple.line3);
		if (triple.tripleBond == null) {
			return;
		}
		JoinPoint firstTripleBond = triple.tripleBond.getJoinPoints().get(0);
		JoinPoint secondTripleBond = triple.tripleBond.getJoinPoints().get(1);
		JoinPoint firstHatchedBond = triple.hatchedBond.getJoinPoints().get(0);
		JoinPoint secondHatchedBond = triple.hatchedBond.getJoinPoints().get(1);
		if ((joinPointsGroupedIntoJunctions.getSizeOfPartition(firstHatchedBond) > 1 || joinPointsGroupedIntoJunctions.getSizeOfPartition(secondHatchedBond) > 1) && joinPointsGroupedIntoJunctions.getSizeOfPartition(firstTripleBond) == 1 && joinPointsGroupedIntoJunctions.getSizeOfPartition(secondTripleBond) == 1) {
			joinPointsGroupedIntoJunctions.remove(firstTripleBond);
			joinPointsGroupedIntoJunctions.remove(secondTripleBond);
			removeJoinable(triple.tripleBond);
		} else {
			joinPointsGroupedIntoJunctions.remove(firstHatchedBond);
			joinPointsGroupedIntoJunctions.remove(secondHatchedBond);
			removeJoinable(triple.hatchedBond);
		}
	}

	private void handleMutuallyExclusiveLineGroupPair(UnionFind<JoinPoint> joinPointsGroupedIntoJunctions, Set<SingleBond> singleBonds, LineGroupPair pair) {
		JoinPoint firstUnsaturatedBond = pair.getUnsaturatedBond().getJoinPoints().get(0);
		JoinPoint secondUnsaturatedBond = pair.getUnsaturatedBond().getJoinPoints().get(1);
		boolean sewn = false;
		for (SingleBond s : pair.getSingleBonds()) {
			try {
				sewn |= joinPointsGroupedIntoJunctions.get(firstUnsaturatedBond).equals(joinPointsGroupedIntoJunctions.get(s.getJoinPoints().get(0)));
				sewn |= joinPointsGroupedIntoJunctions.get(firstUnsaturatedBond).equals(joinPointsGroupedIntoJunctions.get(s.getJoinPoints().get(1)));
				sewn |= joinPointsGroupedIntoJunctions.get(secondUnsaturatedBond).equals(joinPointsGroupedIntoJunctions.get(s.getJoinPoints().get(0)));
				sewn |= joinPointsGroupedIntoJunctions.get(secondUnsaturatedBond).equals(joinPointsGroupedIntoJunctions.get(s.getJoinPoints().get(1)));
			} catch (NullPointerException e) {
				
			}
		}
		if (sewn) {
			Set<JoinPoint> points = joinPointsGroupedIntoJunctions.getObjectsInPartitionOf(firstUnsaturatedBond);
			boolean foundParallel = false;
			for (JoinPoint p1 : points) {
				if (!(p1.getJoinable().getClass().isInstance(pair.getUnsaturatedBond())) && !singleBonds.contains(p1.getJoinable()) && Joinable.areParallel(p1.getJoinable(), pair.getUnsaturatedBond(), new Angle(chemistryBuilder.getParameters().getToleranceForParallelJoinables(), Units.RADIANS))) {
					if (foundParallel) {
						joinPointsGroupedIntoJunctions.explode(points);
						joinPointsGroupedIntoJunctions.remove(firstUnsaturatedBond);
						joinPointsGroupedIntoJunctions.remove(secondUnsaturatedBond);
						removeJoinable(pair.getUnsaturatedBond());
						Set<Joinable> joinables = new HashSet<Joinable>();
						for (JoinPoint p2 : points) {
							if (p2.getJoinable() != pair.getUnsaturatedBond()) {
								joinables.add(p2.getJoinable());
							}
						}
						chemistryBuilder.attemptToJoinListOfJoinables(new ArrayList<Joinable>(joinables), joinPointsGroupedIntoJunctions);
						return;
					} else {
						foundParallel = true;
					}
				}
			}
			points = joinPointsGroupedIntoJunctions.getObjectsInPartitionOf(secondUnsaturatedBond);
			foundParallel = false;
			for (JoinPoint p1 : points) {
				if (!(p1.getJoinable().getClass().isInstance(pair.getUnsaturatedBond())) && !singleBonds.contains(p1.getJoinable()) && Joinable.areParallel(p1.getJoinable(), pair.getUnsaturatedBond(), new Angle(chemistryBuilder.getParameters().getToleranceForParallelJoinables(), Units.RADIANS))) {
					if (foundParallel) {
						joinPointsGroupedIntoJunctions.explode(points);
						joinPointsGroupedIntoJunctions.remove(firstUnsaturatedBond);
						joinPointsGroupedIntoJunctions.remove(secondUnsaturatedBond);
						removeJoinable(pair.getUnsaturatedBond());
						Set<Joinable> joinables = new HashSet<Joinable>();
						for (JoinPoint p2 : points) {
							if (p2.getJoinable() != pair.getUnsaturatedBond()) {
								joinables.add(p2.getJoinable());
							}
						}
						chemistryBuilder.attemptToJoinListOfJoinables(new ArrayList<Joinable>(joinables), joinPointsGroupedIntoJunctions);
						return;
					} else {
						foundParallel = true;
					}
				}
			}
			for (SingleBond s : pair.getSingleBonds()) {
				try {
					joinPointsGroupedIntoJunctions.remove(s.getJoinPoints().get(0));
					joinPointsGroupedIntoJunctions.remove(s.getJoinPoints().get(1));
					removeJoinable(s);
				} catch (NullPointerException e) {
					
				}
			}
			for (SVGLine l : pair.getLines()) {
				chemistryBuilder.getHigherPrimitives().getLineList().remove(l);
			}
		} else {
			for (SingleBond s : pair.getSingleBonds()) {
				try {
					joinPointsGroupedIntoJunctions.remove(s.getJoinPoints().get(0));
					joinPointsGroupedIntoJunctions.remove(s.getJoinPoints().get(1));
					removeJoinable(s);
				} catch (NullPointerException e) {
					
				}
			}
			for (SVGLine l : pair.getLines()) {
				chemistryBuilder.getHigherPrimitives().getLineList().remove(l);
			}
		}
	}
	
	private void undoDamageFromIncorrectMinus(UnionFind<JoinPoint> joinPointsGroupedIntoJunctions, JoinPoint minus) {
		if (minus != null) {
			Set<JoinPoint> points = joinPointsGroupedIntoJunctions.getObjectsInPartitionOf(minus);
			joinPointsGroupedIntoJunctions.explode(points);
			Set<Joinable> joinables = new HashSet<Joinable>();
			for (JoinPoint p : points) {
				if (p.getJoinable() != minus.getJoinable()) {
					joinables.add(p.getJoinable());
				}
			}
			chemistryBuilder.attemptToJoinListOfJoinables(new ArrayList<Joinable>(joinables), joinPointsGroupedIntoJunctions);
		}
	}

	boolean mutuallyExclusive(Joinable joinableI, Joinable joinableJ) {
		for (MutuallyExclusiveShortLineTriple triple : mutuallyExclusiveShortLineTriples) {
			if (joinableI == triple.hatchedBond && joinableJ == triple.minus || joinableJ == triple.hatchedBond && joinableI == triple.minus) {
				return true;
			}
			if (joinableI == triple.hatchedBond && joinableJ == triple.singleBond || joinableJ == triple.hatchedBond && joinableI == triple.singleBond) {
				return true;
			}
			if (joinableI == triple.singleBond && joinableJ == triple.minus || joinableJ == triple.singleBond && joinableI == triple.minus) {
				return true;
			}
		}
		
		for (MutuallyExclusiveShortLinePairTriple triple : mutuallyExclusiveShortLinePairTriples) {
			if (joinableI == triple.hatchedBond && joinableJ == triple.doubleBond || joinableJ == triple.hatchedBond && joinableI == triple.doubleBond) {
				return true;
			}
			if (joinableI == triple.hatchedBond && joinableJ == triple.singleBond1 || joinableJ == triple.hatchedBond && joinableI == triple.singleBond1) {
				return true;
			}
			if (joinableI == triple.singleBond1 && joinableJ == triple.doubleBond || joinableJ == triple.singleBond1 && joinableI == triple.doubleBond) {
				return true;
			}
			if (joinableI == triple.hatchedBond && joinableJ == triple.singleBond2 || joinableJ == triple.hatchedBond && joinableI == triple.singleBond2) {
				return true;
			}
			if (joinableI == triple.doubleBond && joinableJ == triple.singleBond2 || joinableJ == triple.doubleBond && joinableI == triple.singleBond2) {
				return true;
			}
			if (joinableI == triple.singleBond1 && joinableJ == triple.singleBond2 || joinableJ == triple.singleBond1 && joinableI == triple.singleBond2) {
				return true;
			}
		}
		
		/*for (MutuallyExclusiveLinePairPair pair : mutuallyExclusiveLinePairPairs) {
			if (joinableI == pair.singleBond2 && joinableJ == pair.doubleBond || joinableJ == pair.singleBond2 && joinableI == pair.doubleBond) {
				return true;
			}
			if (joinableI == pair.singleBond2 && joinableJ == pair.singleBond1 || joinableJ == pair.singleBond2 && joinableI == pair.singleBond1) {
				return true;
			}
			if (joinableI == pair.singleBond1 && joinableJ == pair.doubleBond || joinableJ == pair.singleBond1 && joinableI == pair.doubleBond) {
				return true;
			}
		}*/
		
		return false;
	}
	
}