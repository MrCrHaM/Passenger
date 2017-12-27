Similarity-based Clustering
---------------------------

A classic [tutorial](http://www.kyb.mpg.de/fileadmin/user_upload/files/publications/attachments/luxburg06_TR_v2_4139%5b1%5d.pdf) on spectral clustering by von Luxburg: this is a great reference, with details of a few versions of the method and pseudocode.

Here is a paper of mine that deals with the semi-supervised setting in which we are looking for a cluster around a certain vertex in a graph. This can be useful if we have detected an anomalous pixel and we want to extend it to an anomalous region: [A Local Spectral Method for Graphs: With Applications to Improving
Graph Partitions and Exploring Data Graphs Locally](http://jmlr.csail.mit.edu/papers/volume13/mahoney12a/mahoney12a.pdf) Notice that the algorithm just runs PageRank on the graph.

Here is some sample code in Matlab to give you an idea. 