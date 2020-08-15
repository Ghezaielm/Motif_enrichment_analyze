# Motif_enrichment_analyze
A concise python class to perform motif enrichment analysis using Python 

## Overview
Let <b>Mm</b> be a PFM for a given motif of length <b>l</b>, <b>S</b> be a query sequence of length <b>L</b> and <b>Ms(i,j)</b> be a PFM of the sub string <b>S[i:j]</b>.<br>

<b><i>The method measure motif enrichment as following:</i></b><br>

For each sub string of <b>S</b> with length <b>l</b>:<br>

a) Compute the dot product for an hypothetical perfect match: Mm.PFM(max occuring NT) = Mp, with size l,l.<br>
b) Compute the dot product for an hypothetical worst match: Mm.PFM(min occuring NT) = Mw, with size l,l.<br>
c) Finally, compute the dot product for the query match: Mm.Ms(i,j) = Mq, with size l,l.<br>

<b>The enrichment score is obtained as following:</b><br>

<b>Enrichment</b> = mean(1-(tr(Mq(i,j)-tr(Mp)-tr(Mw)/(tr/Mp))²/std(1-(tr(Mq(i,j)-tr(Mp)-tr(Mw))/(tr/Mp))

<i>If the query is matching the perfect scenario, then:</i><br>
  <b>Enrichment</b> = mean(1-(-tr(Mw))/(tr(Mp)))²/std(1-(tr(Mw))/(tr/Mp))
 
<i>In a non uniform nucleotide distribution at a single position, tr(Mw) is expected to be very small compared to the perfect match (Mp), so:</i><br>
  <b>Enrichment</b> ~= 1/sigma² with sigma² small, so the results is expected to be above one depending on the standard deviation. 

<i>If the query is matching the worst scenario, then:</i><br>
<b>Enrichment</b> = mean(1-(-tr(Mp))/(tr(Mp)))²/std<br>
So the results will be near to zero.


## Input motifs:
I used the TF2DNA database of human TFs as input for the algorithm. 
