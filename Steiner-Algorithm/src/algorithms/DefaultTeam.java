package algorithms;

import java.awt.Point;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import algorithms.DefaultTeam.Edge;
import algorithms.DefaultTeam.TaggedPoint;

public class DefaultTeam {
	//Sort function : 

	private ArrayList<Edge> sort(ArrayList<Edge> edges) {
		if (edges.size()==1) return edges;

		ArrayList<Edge> left = new ArrayList<Edge>();
		ArrayList<Edge> right = new ArrayList<Edge>();
		int n=edges.size();
		for (int i=0;i<n/2;i++) { left.add(edges.remove(0)); }
		while (edges.size()!=0) { right.add(edges.remove(0)); }
		left = sort(left);
		right = sort(right);

		ArrayList<Edge> result = new ArrayList<Edge>();
		while (left.size()!=0 || right.size()!=0) {
			if (left.size()==0) { result.add(right.remove(0)); continue; }
			if (right.size()==0) { result.add(left.remove(0)); continue; }
			if (left.get(0).distance() < right.get(0).distance()) result.add(left.remove(0));
			else result.add(right.remove(0));
		}
		return result;
	}

	private boolean isLeaf(ArrayList<Edge>list, Point p) {
		ArrayList<Edge> cpt = new ArrayList<Edge>();

		for(Edge tmp : list) {
			if(tmp.p == p || tmp.q == p) {
				cpt.add(tmp);
			}
			if(cpt.size()==2) {
				return false;
			}
		}
		return true;

	}
	/*private boolean isLeaf(ArrayList<Edge>list, Edge e) {
		for(Edge tmp : list) {
			if(tmp==e) {
				continue;
			}

			if(tmp.p==e.p || tmp.q==e.p) {
				for(Edge tmp2 : list) {
					if(tmp2==e) {
						continue;
					}
					if(tmp2.p==e.q || tmp2.q==e.q) {
						return false;
					}
				}
			}
		}
		return true;

	}*/


	public int distance(Point p, Point q) {
		return  (int)(Math.sqrt((p.x - q.x)*(p.x - q.x) + (p.y - q.y)*(p.y - q.y)));

	}



	public class Edge{
		private Point p;
		private Point q;

		Edge(Point p, Point q){
			this.p = p;
			this.q = q;
		}

		public Point getP() {
			return p;
		}

		public Point getQ() {
			return q;
		}

		public double distance() {
			return  (Math.sqrt((p.x - q.x)*(p.x - q.x) + (p.y - q.y)*(p.y - q.y)));

		}
	}


	public class TaggedPoint{
		Point p;
		int tag;

		public TaggedPoint(Point p, int tag) {
			this.p = p;
			this.tag = tag;
		}

		public int getTag() {
			return tag;
		}

		public void setTag(int tag) {
			this.tag = tag;
		}

		public Point getP() {
			return p;
		}
	}

	public boolean isSolution(Edge e, ArrayList<TaggedPoint> l) {
		Point p1 = e.p;
		Point p2 = e.q;
		int x=Integer.MIN_VALUE;
		int y=Integer.MIN_VALUE;
		int cpt = 0;
		for (TaggedPoint tp : l) {
			if(tp.p.equals(p1)) {
				x = tp.tag;
				cpt++;
			}
			if(tp.p.equals(p2)) {
				y = tp.tag;
				cpt++;
			}
			if(cpt == 2) {
				break;
			}
		}
		if(x == y) {
			return false;
		}
		for (TaggedPoint tp : l) {
			if(tp.tag == x) {
				tp.tag = y;
			}
		}
		return true;

	}


	public ArrayList<Tree2D> edgesToTree(ArrayList<Edge> l, Point root) {

		ArrayList<Tree2D> subtree = new ArrayList<Tree2D>();
		ArrayList<Edge> lbis = new ArrayList<Edge>(l);
		//lbis.addAll(l);
		for(Edge e : l) {
			if(e.p == root) {
				lbis.remove(e);
				if(!lbis.isEmpty()) {
					subtree.add(new Tree2D(e.q, edgesToTree(lbis, e.q)));	  
				}else{
					subtree.add(new Tree2D(e.q, new ArrayList<Tree2D>()));
				}
			}
			if(e.q == root) {
				lbis.remove(e);
				if(!lbis.isEmpty()) {
					subtree.add(new Tree2D(e.p, edgesToTree(lbis, e.p)));

				}else {
					subtree.add(new Tree2D(e.p, new ArrayList<Tree2D>()));

				}
			}
		}
		return subtree;

	}


	// Calcul du chemin minimum d'un point a un autre

	public double dist(ArrayList<Point> points,int[][] paths, int i, int j) {
		if(paths[i][j]==i) {
			return Point.distance(points.get(i).x, points.get(i).y, points.get(paths[i][j]).x, points.get(paths[i][j]).y);
		}else {
			double d = Point.distance(points.get(i).x, points.get(i).y, points.get(paths[i][j]).x, points.get(paths[i][j]).y);
			return d + dist(points,paths, paths[i][j], j );
		}

	}

	/**
	 * 
	 * @param points set of all points
	 * @param paths
	 * @param i first point index
	 * @param j second point index
	 * @return list of indexes of points traveled to go from i to j
	 */



	public ArrayList<Integer> pathFinder(ArrayList<Point> points,int[][] paths, int i, int j) {
		ArrayList<Integer> liste = new ArrayList<Integer>();
		while(paths[i][j]!=j) {
			liste.add(paths[i][j]);
			int tmp = paths[i][j];
			i=j;
			j=tmp;
		}
		return liste;

	}

	public int[][] calculShortestPaths(ArrayList<Point> points, int edgeThreshold) {
		int[][] paths=new int[points.size()][points.size()];
		double[][] distances=new double[points.size()][points.size()];
		for(int i=0;i<distances.length;i++) {
			for(int j=0;j<distances.length;j++) {
				double d = (Point.distance(points.get(i).x, points.get(i).y, points.get(j).x, points.get(j).y));
				if(d>= edgeThreshold) {
					distances[i][j] = Double.MAX_VALUE;
					paths[i][j] = j;
				}else{
					distances[i][j] = d;
					paths[i][j] = j;
				}
			}
		}
		for (int k=0;k<paths.length;k++) {
			for (int i=0;i<paths.length;i++) {
				for(int j=0;j<paths.length;j++) {
					double d = distances[i][k]+distances[k][j];
					if(d < distances[i][j] && d > 0) {
						distances[i][j] = d;
						paths[i][j] = paths[i][k];
					}
				}
			}
		}


		return paths;
	}
	// __________________________________________ //


	public Tree2D calculSteiner(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		int[][] paths = calculShortestPaths(points, edgeThreshold);


		int[] indexHitPoint = new int[hitPoints.size()];

		for (int i = 0; i<hitPoints.size(); i++) {
			indexHitPoint[i]= points.indexOf(hitPoints.get(i));
		}

		//First Kruskal

		ArrayList<Point> points2 = new ArrayList<Point>(hitPoints);
		ArrayList<Edge> edges = new ArrayList<Edge>();
		for(int i = 0; i < hitPoints.size(); i++)
			for(int j = i+1; j < hitPoints.size(); j++)
				edges.add(new Edge(hitPoints.get(i), hitPoints.get(j)));

		edges = sort(edges);

		ArrayList<Edge> solution = new ArrayList<Edge>();
		ArrayList<TaggedPoint> TaggedP = new ArrayList<TaggedPoint>();

		for(int i = 0; i < hitPoints.size(); i++) {
			TaggedP.add(new TaggedPoint(hitPoints.get(i), i));
		}

		for(Edge a :edges) {
			if(points2.isEmpty()) {
				break;
			}
			if( isSolution(a, TaggedP) ) {
				points2.remove(a.p);
				solution.add(a);
			}
		}

		// End

		//Replacement of the edge by a shorter path


		for(Edge e : solution) {
			Point tmp1 = e.p;
			Point tmp2 = e.q;
			int i1 = points.indexOf(tmp1);
			int i2 = points.indexOf(tmp2);
			ArrayList<Integer> l = pathFinder(points, paths, i1, i2); 
			for(int i=0; i<l.size(); i++) {
				if(!(hitPoints.contains(points.get(l.get(i))))){
					hitPoints.add(points.get(l.get(i)));
				}
			}  
		}

		points2 = new ArrayList<Point>(hitPoints);
		edges = new ArrayList<Edge>();
		for(int i = 0; i < hitPoints.size(); i++)
			for(int j = i+1; j < hitPoints.size(); j++)
				edges.add(new Edge(hitPoints.get(i), hitPoints.get(j)));


		// Second Kruskal




		edges = sort(edges);

		solution = new ArrayList<Edge>();
		TaggedP = new ArrayList<TaggedPoint>();

		for(int i = 0; i < hitPoints.size(); i++) {
			TaggedP.add(new TaggedPoint(hitPoints.get(i), i));
		}
		int total = 0;
		for(Edge e :edges) {
			if(points2.isEmpty() ) {
				break;
			}
			if( isSolution(e, TaggedP) ) {
				total+= dist(points, paths, points.indexOf(e.p), points.indexOf(e.q));
				points2.remove(e.p);
				solution.add(e);
			}
		}

		// Solution
		double cost = 0;
		for(Edge e: solution) {
			cost+= e.distance();
		}
		System.out.println(cost);
		return new Tree2D(solution.get(0).p, edgesToTree(solution, solution.get(0).p));


	}

	//__________________________Steiner width budget______________________________//

	public Tree2D calculSteinerBudget(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
		int[][] paths = calculShortestPaths(points, edgeThreshold);
		ArrayList<Point> savedHP = new ArrayList<Point>(hitPoints);

		int[] indexHitPoint = new int[hitPoints.size()];

		for (int i = 0; i<hitPoints.size(); i++) {
			indexHitPoint[i]= points.indexOf(hitPoints.get(i));
		}

		//First Kruskal


		ArrayList<Point> points2 = new ArrayList<Point>(hitPoints);
		ArrayList<Edge> edges = new ArrayList<Edge>();	
		for(int i = 0; i < hitPoints.size(); i++)
			for(int j = i+1; j < hitPoints.size(); j++)
				edges.add(new Edge(hitPoints.get(i), hitPoints.get(j)));


		edges = sort(edges);


		ArrayList<Edge> solution = new ArrayList<Edge>();
		ArrayList<TaggedPoint> TaggedP = new ArrayList<TaggedPoint>();

		for(int i = 0; i < hitPoints.size(); i++) {
			TaggedP.add(new TaggedPoint(hitPoints.get(i), i));
		}

		for(Edge a :edges) {
			if(points2.isEmpty()) {
				break;
			}
			if( isSolution(a, TaggedP) ) {
				points2.remove(a.p);
				solution.add(a);
			}
		}



		// End

		//Replacement of the edge by a shorter path


		for(Edge e : solution) {
			Point tmp1 = e.p;
			Point tmp2 = e.q;
			int i1 = points.indexOf(tmp1);
			int i2 = points.indexOf(tmp2);
			ArrayList<Integer> l = pathFinder(points, paths, i1, i2); 
			for(int i=0; i<l.size(); i++) {
				if(!(hitPoints.contains(points.get(l.get(i))))){
					hitPoints.add(points.get(l.get(i)));
				}
			}  
		}

		points2 = new ArrayList<Point>(hitPoints);
		edges = new ArrayList<Edge>();
		for(int i = 0; i < hitPoints.size(); i++)
			for(int j = i+1; j < hitPoints.size(); j++)
				edges.add(new Edge(hitPoints.get(i), hitPoints.get(j)));


		// Second Kruskal




		edges = sort(edges);

		solution = new ArrayList<Edge>();
		TaggedP = new ArrayList<TaggedPoint>();

		for(int i = 0; i < hitPoints.size(); i++) {
			TaggedP.add(new TaggedPoint(hitPoints.get(i), i));
		}
		for(Edge e :edges) {
			if(points2.isEmpty() ) {
				break;
			}
			if( isSolution(e, TaggedP) ) {
				//cost+= dist(points, paths, points.indexOf(e.p), points.indexOf(e.q));
				points2.remove(e.p);
				solution.add(e);
			}
		}
		//End kruskal_______________________________________________//


		double cost = 0;
		for(Edge e: solution) {
			cost+= e.distance();
		}

		Edge root = solution.get(0);
		ArrayList<Edge> removedE = new ArrayList<Edge>();
		ArrayList<Edge> copySol = new ArrayList<Edge>(solution);
		System.out.println(cost);

		while(cost > 1664) {
			Point MostFar = null;
			double distance = 0;

			for (Point p : hitPoints) {
				if(distance(root.p,p)>distance) {
					distance = dist(points, paths, points.indexOf(root.p), points.indexOf(p));
					MostFar = p;
				}
				if(distance(root.q,p)>distance) {
					distance = dist(points, paths, points.indexOf(root.q), points.indexOf(p));
					MostFar = p;
				}
			}
			Edge toRemove = null;

			for(Edge e : solution) {
				if((isLeaf(solution,e.p) && e.p == MostFar) || isLeaf(solution,e.q) &&  e.q == MostFar) {
					cost -= e.distance();
					toRemove = e;
				}
			}
			solution.remove(toRemove);
			copySol.remove(toRemove);
			hitPoints.remove(MostFar);
			boolean aSuppr = true;

			while(aSuppr) { //delete useless edges
				boolean test = false;
				Edge tmpe = null;
				for(Edge e : copySol) {
					if(isLeaf(copySol,e.p)|| isLeaf(copySol,e.q)) {
						if(!savedHP.contains(e.p) || !savedHP.contains(e.q)  ) {
							test = true;
							tmpe = e;
							break;
						}
					}
				}
				if(test) {
					removedE.add(tmpe);
					copySol.remove(tmpe);
					cost -= tmpe.distance();
				}else {
					break;
				}
			}

		}	

		double costsol2 = 0;
		for(Edge e: copySol) {
			costsol2+=e.distance();
		}

		System.out.println("\n"+costsol2);
		
		return new Tree2D(copySol.get(0).p, edgesToTree(copySol, copySol.get(0).p));

	}


}





