// This file is required in order for any other classes to work. Some Vector methods work with the
// other Sylvester classes and are useless unless they are included. Other classes such as Line and
// Plane will not function at all without Vector being loaded first.

var Sylvester = {
  precision: 1e-6
};
function Vector() {}
Vector.prototype = {

  // Returns element i of the vector
  e: function e(i) {
    return (i < 1 || i > this.elements.length) ? null : this.elements[i-1];
  },
  
  // Returns element i of the sorted vector
  eSorted: function eSorted(i) {
	if (!this.sorted) {this.saveSorted();}
	return (i < 1 || i > this.sortedElements.length) ? null : this.sortedElements[i-1];
  },
  
  // Edits element i of the vector
  set: function set(i, val) {
    if (i < 1 || i > this.elements.length) { return null; };
	this.elements[i-1] = val;
	this.sorted=false;
  },

  // Returns the number of elements the vector has
  dimensions: function dimensions() {
    return this.elements.length;
  },

  // Returns the modulus ('length') of the vector
  modulus: function modulus() {
    return Math.sqrt(this.dot(this));
  },

  // Returns true iff the vector is equal to the argument
  eql: function eql(vector) {
    var n = this.elements.length;
    var V = vector.elements || vector;
    if (n != V.length) { return false; }
	while (n--) {
      if (Math.abs(this.elements[n] - V[n]) > Sylvester.precision) { return false; }
    }
    return true;
  },

  // Returns a copy of the vector
  dup: function dup() {
    return Vector.create(this.elements);
  },

  // Maps the vector to another vector according to the given function
  map: function map(fn) {
    var elements = [];
    this.each(function(x, i) {
      elements.push(fn(x, i));
    });
    return Vector.create(elements);
  },
  
  // Calls the iterator for each element of the vector in turn
  each: function each(fn) {
    var n = this.elements.length, k = n, i;
	for (var i = 0; i < n; i++) {
      fn(this.elements[i], i+1);
    }
	this.sorted=false;
  },

  // Returns a new vector created by normalizing the receiver
  toUnitVector: function toUnitVector() {
    var r = this.modulus();
    if (r === 0) { return this.dup(); }
    return this.map(function(x) { return x/r; });
  },

  // Returns the angle between the vector and the argument (also a vector)
  angleFrom: function angleFrom(vector) {
    var V = vector.elements || vector;
    var n = this.elements.length, k = n, i;
    if (n != V.length) { return null; }
    var dot = 0, mod1 = 0, mod2 = 0;
    // Work things out in parallel to save time
    this.each(function(x, i) {
      dot += x * V[i-1];
      mod1 += x * x;
      mod2 += V[i-1] * V[i-1];
    });
    mod1 = Math.sqrt(mod1); mod2 = Math.sqrt(mod2);
    if (mod1*mod2 === 0) { return null; }
    var theta = dot / (mod1*mod2);
    if (theta < -1) { theta = -1; }
    if (theta > 1) { theta = 1; }
    return Math.acos(theta);
  },

  // Returns true iff the vector is parallel to the argument
  isParallelTo: function isParallelTo(vector) {
    var angle = this.angleFrom(vector);
    return (angle === null) ? null : (angle <= Sylvester.precision);
  },

  // Returns true iff the vector is antiparallel to the argument
  isAntiparallelTo: function isAntiparallelTo(vector) {
    var angle = this.angleFrom(vector);
    return (angle === null) ? null : (Math.abs(angle - Math.PI) <= Sylvester.precision);
  },

  // Returns true iff the vector is perpendicular to the argument
  isPerpendicularTo: function isPerpendicularTo(vector) {
    var dot = this.dot(vector);
    return (dot === null) ? null : (Math.abs(dot) <= Sylvester.precision);
  },

  // Returns the result of adding the argument to the vector
  add: function add(vector) {
    var V = vector.elements || vector;
    if (this.elements.length != V.length) { return null; }
    return this.map(function(x, i) { return x + V[i-1]; });
  },

  // Returns the result of subtracting the argument from the vector
  subtract: function subtract(vector) {
    var V = vector.elements || vector;
    if (this.elements.length != V.length) { return null; }
    return this.map(function(x, i) { return x - V[i-1]; });
  },

  // Returns the result of multiplying the elements of the vector by the argument
  multiply: function multiply(k) {
    return this.map(function(x) { return x*k; });
  },

  x: function x(k) { return this.multiply(k); },

  // Returns the scalar product of the vector with the argument
  // Both vectors must have equal dimensionality
  dot: function dot(vector) {
    var V = vector.elements || vector;
    var i, product = 0, n = this.elements.length;
    if (n != V.length) { return null; }
    do { product += this.elements[n-1] * V[n-1]; } while (--n);
    return product;
  },

  // Returns the vector product of the vector with the argument
  // Both vectors must have dimensionality 3
  cross: function cross(vector) {
    var B = vector.elements || vector;
    if (this.elements.length != 3 || B.length != 3) { return null; }
    var A = this.elements;
    return Vector.create([
      (A[1] * B[2]) - (A[2] * B[1]),
      (A[2] * B[0]) - (A[0] * B[2]),
      (A[0] * B[1]) - (A[1] * B[0])
    ]);
  },
  
  // Saves a sorted equivalent of the vector
  saveSorted: function saveSorted() {
	if (!this.sorted) {
		this.sortedElements = this.elements.slice();
		this.sortedElements.sort(function callbackFunc(a, b){
	return a - b;
});
		this.sorted = true;
	}
  },
  
  // Returns the sorted array
  sort: function sort() {
	if (!this.sorted) {
		this.saveSorted();
	}
	return this.sortedElements;
  },

  // Returns the (absolute) largest element of the vector
  maxAbs: function maxAbs() {
    var m = 0, i = this.elements.length;
    while (i--) {
      if (Math.abs(this.elements[i]) > Math.abs(m)) { m = this.elements[i]; }
    }
  },
  
  max: function max() {
	if (!this.sorted) {
		this.saveSorted();
	}
	return this.sortedElements[this.sortedElements.length];
    // var m = 0, n = this.elements.length, k = n, i;
    // do { i = k - n;
      // if (this.elements[i] > m) { m = this.elements[i]; }
    // } while (--n);
    // return m;
  },
  
  max2: function max2() {
    var m = 0, i = this.elements.length;
    while (i--) {
      if (this.elements[i] > m) { m = this.elements[i]; }
    }
    return m;
  },
  
  // Returns the smallest element of the vector
  min: function min() {
	if (!this.sorted) {
		this.saveSorted();
	}
	return this.sortedElements[0];
    // var m = 0, n = this.elements.length, k = n, i;
    // do { i = k - n;
      // if (this.elements[i] < m) { m = this.elements[i]; }
    // } while (--n);
    // return m;
  },
  min2: function min2() {
    var m = 0, n = this.elements.length, k = n, i;
    do { i = k - n;
      if (this.elements[i] < m) { m = this.elements[i]; }
    } while (--n);
    return m;
  },

  // Returns the index of the first match found
  indexOf: function indexOf(x) {
    var n = this.elements.length, i;
	for (i=0; i<n; i++) {
		if (this.elements[i] == x) { return i+1; }
	}
    return null;
  },

  // Returns a diagonal matrix with the vector's elements as its diagonal elements
  toDiagonalMatrix: function toDiagonalMatrix() {
    return Matrix.Diagonal(this.elements);
  },

  // Returns the result of rounding the elements of the vector
  round: function round() {
    return this.map(function(x) { return Math.round(x); });
  },

  // Returns a copy of the vector with elements set to the given value if they
  // differ from it by less than Sylvester.precision
  snapTo: function snapTo(x) {
    return this.map(function(y) {
      return (Math.abs(y - x) <= Sylvester.precision) ? x : y;
    });
  },

  // Returns the vector's distance from the argument, when considered as a point in space
  distanceFrom: function distanceFrom(obj) {
    if (obj.anchor) { return obj.distanceFrom(this); }
    var V = obj.elements || obj;
    if (V.length != this.elements.length) { return null; }
    var sum = 0, part;
    this.each(function(x, i) {
      part = x - V[i-1];
      sum += part * part;
    });
    return Math.sqrt(sum);
  },

  // Returns true if the vector is point on the given line
  liesOn: function liesOn(line) {
    return line.contains(this);
  },

  // Return true iff the vector is a point in the given plane
  liesIn: function liesIn(plane) {
    return plane.contains(this);
  },

  // Rotates the vector about the given object. The object should be a 
  // point if the vector is 2D, and a line if it is 3D. Be careful with line directions!
  rotate: function rotate(t, obj) {
	var V, R = null, x, y, z;
    if (t.determinant) { R = t.elements; }
    switch (this.elements.length) {
      case 2:
        V = obj.elements || obj;
        if (V.length != 2) { return null; }
        if (!R) { R = Matrix.Rotation(t).elements; }
        x = this.elements[0] - V[0];
        y = this.elements[1] - V[1];
        return Vector.create([
          V[0] + R[0][0] * x + R[0][1] * y,
          V[1] + R[1][0] * x + R[1][1] * y
        ]);
        break;
      case 3:
        if (!obj.direction) { return null; }
        var C = obj.pointClosestTo(this).elements;
        if (!R) { R = Matrix.Rotation(t, obj.direction).elements; }
        x = this.elements[0] - C[0];
        y = this.elements[1] - C[1];
        z = this.elements[2] - C[2];
        return Vector.create([
          C[0] + R[0][0] * x + R[0][1] * y + R[0][2] * z,
          C[1] + R[1][0] * x + R[1][1] * y + R[1][2] * z,
          C[2] + R[2][0] * x + R[2][1] * y + R[2][2] * z
        ]);
        break;
      default:
        return null;
    }
  },

  // Returns the result of reflecting the point in the given point, line or plane
  reflectionIn: function reflectionIn(obj) {
    if (obj.anchor) {
      // obj is a plane or line
      var P = this.elements.slice();
      var C = obj.pointClosestTo(P).elements;
      return Vector.create([C[0] + (C[0] - P[0]), C[1] + (C[1] - P[1]), C[2] + (C[2] - (P[2] || 0))]);
    } else {
      // obj is a point
      var Q = obj.elements || obj;
      if (this.elements.length != Q.length) { return null; }
      return this.map(function(x, i) { return Q[i-1] + (Q[i-1] - x); });
    }
  },

  // Utility to make sure vectors are 3D. If they are 2D, a zero z-component is added
  to3D: function to3D() {
    var V = this.dup();
    switch (V.elements.length) {
      case 3: break;
      case 2: V.elements.push(0); break;
      default: return null;
    }
    return V;
  },

  // Returns a string representation of the vector
  inspect: function inspect() {
    return '[' + this.elements.join(', ') + ']';
  },

  // Set vector's elements from an array
  setElements: function setElements(els) {
    this.elements = (els.elements || els).slice();
	this.sorted=false;
    return this;
  }
};
  
// Constructor function
Vector.create = function VectorCreate(elements) {
  var V = new Vector();
  
  return V.setElements(elements);
};

// i, j, k unit vectors
Vector.i = Vector.create([1,0,0]);
Vector.j = Vector.create([0,1,0]);
Vector.k = Vector.create([0,0,1]);

// Random vector of size n
Vector.Random = function VectorRandom(n) {
  var elements = [];
  while (n--) { elements.push(Math.random()); }
  return Vector.create(elements);
};

// Vector filled with zeros
Vector.Zero = function VectorZero(n) {
  var elements = [];
  while (n--) { elements.push(0); }
  return Vector.create(elements);
};