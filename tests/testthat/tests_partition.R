library(testthat)


dag = create_ontology_DAG(
	c("a - b",
	  "a - c",
	  "a - d",
	  "b - e",
	  "b - f",
	  "c - g", 
	  "d - h",
	  "g - i",
	  "g - j",
	  "h - k",
	  "h - l",
	  "i - m",
	  "j - n",
	  "k - o",
	  "l - p")
)

test_that("test partition_by_level", {
	expect_equal(
		partition_by_level(dag, level = 0),
		rep("a", 16)
	)
	expect_equal(
		partition_by_level(dag, level = 1),
		c(NA, "b", "c", "d", "b", "b", "c", "d", "c", "c", "d", "d", "c", "c", "d", "d")
	)
	expect_equal(
		partition_by_level(dag, level = 2),
		c(NA, NA, NA, NA, "e", "f", "g", "h", "g", "g", "h", "h", "g", "g", "h", "h")
	)
	
})

test_that("test partition_by_size", {
	expect_equal(
		partition_by_size(dag, size = 3),
		c(NA, "b", NA, NA, "b", "b", "g", "h", "g", "g", "h", "h", "g", "g", "h", "h")
	)

	expect_equal(
		partition_by_size(dag, size = 6),
		c(NA, "b", "c", "d", "b", "b", "c", "d", "c", "c", "d", "d", "c", "c", "d", "d")
	)
})
