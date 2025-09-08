#
# Only gives the nearest neighbor, not the k-nearest neighbors.
# If you give as query the same data that makes up the tree
# all nearest neighbor distances should be 0. They are not,so there is a mistake.
# There are 3 different routines to insert points to the kd tree
# insertRecursive
# insertRecursive_old
# and bulkInsert
# There does not seem to be a huge speed difference. I think that R is just too
# slow for this.



# Load the R6 library if not already loaded
if (!requireNamespace("R6", quietly = TRUE)) {
  install.packages("R6")
}
library(R6)

# Define the Point class
Point <- R6Class("Point",
                 public = list(
                   u = NULL,
                   v = NULL,
                   initialize = function(u, v) {
                     self$u <- u
                     self$v <- v
                   }
                 ))

# Define the Node class, inheriting from Point
Node <- R6Class(
  "Node",
  inherit = Point,
  public = list(
    left = NULL,
    right = NULL,
    splitDim = "u",
    initialize = function(u = NA,
                          v = NA,
                          left = NULL,
                          right = NULL,
                          splitDim = "u") {
      super$initialize(u = u, v = v)
      self$left <- left
      self$right <- right
      self$splitDim <- splitDim
    },


    insertRecursive = function(point) {

      if (is.na(self$u) && is.na(self$v)) {
        # The tree is empty, insert the point at the root
        self$u <- point$u
        self$v <- point$v
        return()
      }

      current <- self
      nextSplitDim <- self$splitDim

      # Traverse the tree iteratively
      while (!is.null(current)) {
        if ((nextSplitDim == "u" && point$u < current$u) || (nextSplitDim == "v" && point$v < current$v)) {
          if (is.null(current$left)) {
            # Create a new left node and insert the point
            current$left <- Node$new(u = point$u, v = point$v, splitDim = nextSplitDim)
            break
          } else {
            # Move to the left subtree
            current <- current$left
          }
        } else {
          if (is.null(current$right)) {
            # Create a new right node and insert the point
            current$right <- Node$new(u = point$u, v = point$v, splitDim = nextSplitDim)
            break
          } else {
            # Move to the right subtree
            current <- current$right
          }
        }

        # Switch the split dimension for the next level
        nextSplitDim <- if (nextSplitDim == "u") "v" else "u"
      }
    },


    bulkInsert_recursive = function(points, splitDim = "u") {
      if (nrow(points) > 0) {
        self$splitDim <- splitDim

        if (nrow(points) == 1) {
          self$u <- points[1, "u"]
          self$v <- points[1, "v"]
        } else {
          currentDim <- ifelse(splitDim == "u", "v", "u")
          medianIdx <- nrow(points) %/% 2

          sortedIndices <- order(points[, currentDim])
          self$u <- points[sortedIndices[medianIdx], "u"]
          self$v <- points[sortedIndices[medianIdx], "v"]

          leftIndices <- sortedIndices[1:(medianIdx - 1)]
          rightIndices <- sortedIndices[(medianIdx + 1):nrow(points)]

          if (length(leftIndices) > 0) {
            self$left <- Node$new()
            self$left$bulkInsert(points[leftIndices, ], ifelse(splitDim == "u", "v", "u"))
          }
          if (length(rightIndices) > 0) {
            self$right <- Node$new()
            self$right$bulkInsert(points[rightIndices, ], ifelse(splitDim == "u", "v", "u"))
          }
        }
      }
    },
    bulkInsert = function(points, splitDim = "u") {
      if (nrow(points) == 0) {
        return(invisible(self))
      }

      self$splitDim <- splitDim

      stack <- list(node = self, points = points)

      while (length(stack) > 0) {
        current <- stack[[length(stack)]]
        stack <- stack[1:(length(stack) - 1)]

        currentDim <- ifelse(current$splitDim == "u", "v", "u")
        medianIdx <- nrow(current$points) %/% 2

        sortedIndices <- order(current$points[, currentDim])
        current$node$u <- current$points[sortedIndices[medianIdx], "u"]
        current$node$v <- current$points[sortedIndices[medianIdx], "v"]

        leftIndices <- sortedIndices[1:(medianIdx - 1)]
        rightIndices <- sortedIndices[(medianIdx + 1):nrow(current$points)]

        if (length(leftIndices) > 0) {
          current$node$left <- Node$new()
          stack <- c(stack, list(node = current$node$left, points = current$points[leftIndices, ]))
        }
        if (length(rightIndices) > 0) {
          current$node$right <- Node$new()
          stack <- c(stack, list(node = current$node$right, points = current$points[rightIndices, ]))
        }
      }

      invisible(self)
    },


    # Recursive function to insert a point into the KD tree below or at that node.
    insertRecursive_old = function(point) {
      # No data on that node yet.
      if (is.na(self$u) || is.na(self$v)) {
        self$u <- point$u
        self$v <- point$v
      } else {
        #alternating with every layer of the tree, use U or V coordinate to decide
        #if the current point goes to the left (smaller in that corodinate) or right
        #(larger in that coordinate) of the current node
        nextSplitDim <- if (self$splitDim == "u")
          "v"
        else
          "u"
        if ((self$splitDim == "u" &&
             point$u < self$u) || (self$splitDim == "v" && point$v < self$v)) {
          #point will go to the left of the current node
          if (is.null(self$left)) {
            self$left <-
              Node$new(u = point$u,
                       v = point$v,
                       splitDim = nextSplitDim)
          } else {
            self$left$insertRecursive(point)
          }
        } else {
          #point will go to the right of the current node
          if (is.null(self$right)) {
            self$right <-
              Node$new(u = point$u,
                       v = point$v,
                       splitDim = nextSplitDim)
          } else {
            self$right$insertRecursive(point)
          }
        }
      }
    },
    # Recursive function to find the nearest neighbor in the KD tree
    kNearestNeighborsRecursive = function(queryPoint, currentBest = NULL) {
      if (!is.null(self)) {
        # Calculate the distance between the current node and the query point
        dist <-
          sqrt((self$u - queryPoint$u) ^ 2 + (self$v - queryPoint$v) ^ 2)

        # If currentBest is NULL or the current node is closer than currentBest, update currentBest
        if (is.null(currentBest) ||
            dist < currentBest$NodeDist) {
          currentBest <- list(NodeDist = dist, Point = self)
        }

        # Determine which side of the splitting plane to traverse
        nextSplitDim <-
          if (self$splitDim == "u")
            "v"
        else
          "u"
        if ((self$splitDim == "u" &&
             queryPoint$u < self$u) ||
            (self$splitDim == "v" && queryPoint$v < self$v)) {
          # Traverse left if the query point is on the left
          if (!is.null(self$left)) {
            currentBest <-
              self$left$kNearestNeighborsRecursive(queryPoint, currentBest)
          }

          # Check if there could be nearer points on the other side of the splitting plane
          if ((self$splitDim == "u" &&
               abs(queryPoint$u - self$u) < currentBest$NodeDist) ||
              (self$splitDim == "v" &&
               abs(queryPoint$v - self$v) < currentBest$NodeDist)) {
            if (!is.null(self$right)) {
              currentBest <-
                self$right$kNearestNeighborsRecursive(queryPoint, currentBest)
            }
          }
        } else {
          # Traverse right if the query point is on the right
          if (!is.null(self$right)) {
            currentBest <-
              self$right$kNearestNeighborsRecursive(queryPoint, currentBest)
          }

          # Check if there could be nearer points on the other side of the splitting plane
          if ((self$splitDim == "u" &&
               abs(queryPoint$u - self$u) < currentBest$NodeDist) ||
              (self$splitDim == "v" &&
               abs(queryPoint$v - self$v) < currentBest$NodeDist)) {
            if (!is.null(self$left)) {
              currentBest <-
                self$left$kNearestNeighborsRecursive(queryPoint, currentBest)
            }
          }
        }
      }
      currentBest
    }
  )
)

# Define the KDTree class
KDTree <- R6Class(
  "KDTree",
  public = list(
    root = NULL,
    initialize = function() {
      self$root <- Node$new()
    },
    bulkInsert = function(points) {
      self$root$bulkInsert(points)
      invisible(self)
    },
    # Function to insert a point into the KD tree
    insert = function(point) {
      self$root$insertRecursive(point)
      invisible(self)
    },
    insert_old = function(point) {
      self$root$insertRecursive_old(point)
      invisible(self)
    },
    # Function to find the nearest neighbor in the KD tree
    kNearestNeighbors = function(queryPoint) {
      self$root$kNearestNeighborsRecursive(queryPoint)
    },
    # Function to print the tree structure
    printTree = function() {
      self$root$printNode(0)
    }
  )
)

# Add a printNode method to the Node class
Node$set("public", "printNode", function(indent = 0) {
  cat(
    paste0(rep("  ", indent)),
    paste0(
      "SplitDim:",
      self$splitDim,
      ", Point:",
      self$u,
      ",",
      self$v,
      "\n"
    )
  )
  if (!is.null(self$left)) {
    self$left$printNode(indent + 1)
  }
  if (!is.null(self$right)) {
    self$right$printNode(indent + 1)
  }
})









###################################################################
# Create a KD-tree
kdtree <- KDTree$new()

# Insert a matrix of points
kdtree <- KDTree$new()
points <- data.frame(u = sample(seq(1,10),size = 10,replace = T), v = sample(seq(1,10),size = 10,replace = T))
kdtree$bulkInsert(points)


order(points[, self$splitDim])

nearestneighbor_flo <- function(data, query=1){

  n <- 10
  data <- cbind(runif(n),runif(n))
  points <- data.frame(u=data[,1],v=data[,2])
  kdtree <- KDTree$new()
  tic()
  for(i in 1:nrow(data)){
    kdtree$insert_old(Point$new(data[i,1],data[i,2]))
  }
  toc()
  kdtree <- KDTree$new()
  tic()
 kdtree$bulkInsert(points)
  toc()

  kdtree$printTree()
}

# Example usage:
# Create a KDTree instance
kdtree <- KDTree$new()

# Insert some points into the k-d tree
points <-
  list(Point$new(2, 3),
       Point$new(5, 6),
       Point$new(1, 9),
       Point$new(3, 4))
for (point in points) {
  kdtree$insert(point)
}

# Query k-nearest neighbors for a new point
queryPoint <- Point$new(0.2, 0.5)
neighbors <- rep(0,10)
for(i in 1:10){
  neighbors[i] <- kdtree$kNearestNeighbors( Point$new(data[i,1], data[i,2]))$NodeDist
}
kdtree$kNearestNeighbors(Point$new(data[5,1], data[5,2]))

c(data[5,1], data[5,2])

# Display the result
cat("K-Nearest Neighbors:\n")
for (neighbor in neighbors) {
  cat(
    "Distance:",
    neighbor$NodeDist,
    "Point:",
    neighbor$Point$u,
    neighbor$Point$v,
    "\n"
  )
}


# Display the result
cat("K-Nearest Neighbors:\n")
for (i in seq_along(neighbors)) {
  cat(
    "Distance:",
    neighbors[[i]]$NodeDist,
    "Point:",
    neighbors[[i]]$Point$u,
    neighbors[[i]]$Point$v,
    "\n"
  )
}
