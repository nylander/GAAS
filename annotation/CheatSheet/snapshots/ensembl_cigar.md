
<!DOCTYPE HTML PUBLIC "-//IETF//DTD HTML//EN">
<html>

Copy from <a href="https://web.archive.org/web/20050316014923/http://www.ensembl.org:80/Docs/wiki/html/EnsemblDocs/CigarFormat.html">web.archive.org</a>.

<table border="0" cellpadding="0" cellspacing="0" width="100%">
    <tr>
        <td valign="top" width="20" bgcolor="#e2e2ff">&nbsp;&nbsp;&nbsp;</td>
        <td><table border="0" cellpadding="0" cellspacing="0" width="100%" bgcolor="#e2e2ff">
            <tr>
				<td colspan="3">&nbsp;</td>
			</tr>
			<tr>
                <td><a href="/web/20050210133134/http://www.ensembl.org:80/cgi-bin/wiki/search.pl/EnsemblDocs?text=CigarFormat&amp;options=topic"><h2>CigarFormat</h2></a>
                </td>
				<td width="100%">&nbsp;</td>				
                <td align="right" nowrap><h2>EnsemblDocs <font color="#BC384D">W</font><font color="#E4A327">i</font><font color="#008000">k</font><font color="#000080">i</font>&nbsp;&nbsp;&nbsp;</h2></td>
            </tr>
			</table>
        <table border="0" cellpadding="8" cellspacing="0" width="100%" bordercolor="#FFFFFF">
            <tr>
                <td colspan="2" width="100%"><h3>Description</h3>

<p> 

<a href="/web/20050210133134/http://www.ensembl.org:80/Docs/wiki/html/EnsemblDocs/CigarFormat.html">CigarFormat</a> is one of the output formats which can be generated

by <a href="/web/20050210133134/http://www.ensembl.org:80/Docs/wiki/html/EnsemblDocs/Exonerate.html">Exonerate</a>

<p> 

It is also used in the feature tables in the Ensembl database,

but in an altered form.

<p> 

It is designed to contain the minimal information

necessary for the reconstruction of an alignment.

One alignment is described per line, to allow easy manipulation

with UNIX tools.

<p> 

<b>Cigar</b> is an acronym for <b>C</b>oncise <b>I</b>diosyncratic <b>G</b>apped <b>A</b>lignment <b>R</b>eport.

(There is also the related Sugar - Simple Ungapped Alignment Report,

and Vulgar - Verbose Ugly Labelled Gapped Alignment Report :)

<p> 

<hr>

<p> 

<h3>Format</h3>

<p> 

Cigar format looks like this:

<b>

<pre>
cigar: hs989235.cds 5 468 + hsnfg9.embl 25689 27450 + 1916 M 13 I 1 M 35 I 1 M 4 I 1 M 13 D 1 M 4 I 1 M 115 D 404 M 37 D 1 M 164 I 1 M 12 D 898 M 16 I 1 M 12 I 1 M 21 D 1 M 10
</pre>

</b>

<p> 

The fields are as follows:

<ol>
<li>  query identifier

<li>  query start position

<li>  query stop position

<li>  query strand

<li>  target identifier

<li>  target start position

<li>  target stop position

<li>  target strand

<li>  score

</ol>
<p> 

The remaining fields are in pairs, describing the edit path

throught the alignment.  These contain a M,I,D or N

corresponding to a Match, Insert, Delete or iNtron, followed

by the length.

<p> 

<hr>

<p> 

<h3>Example</h3>

<p> 

Below is an the alignment corresponding to the cigar line show above:

<p> 

<pre>

C4 Alignment display:
  Model: est2genome
  Raw score: 1916
  Aligned positions 5-&gt;468 of query
  Aligned positions 25689-&gt;27450 of target

Query: hs989235.cds
Target: hsnfg9.embl

     6 : AAGCTCANCTTGGACCACCGACTCTCGANTGNNTCGCCGCGGGAGCCGGNTGGANAACCT :    64
         ||||||| ||||| |||||||||||||  ||   ||||||||||||||| |||| |||||
 25690 : AAGCTCATCTTGG-CCACCGACTCTCGCTTGCGCCGCCGCGGGAGCCGG-TGGA-AACCT : 25745

    65 : GAGCGGGA-CTGGNAGAAGGAGCAGAGGGAGGCAGCACCCGGCGTGACGGNAGTGTGTGG :   123
         |||||||| |||| |||||||||||||||||||||||||||||||||||| |||||||||
 25746 : GAGCGGGAGCTGG-AGAAGGAGCAGAGGGAGGCAGCACCCGGCGTGACGGGAGTGTGTGG : 25804

   124 : GGCACTCAGGCCTTCCGCAGTGTCATCTGCCACACGGAAGGCACGGCCACGGGCAGGGGG :   183
         ||||||||||||||||||||||||||||||||||||||||||||||||||||||  ||||
 25805 : GGCACTCAGGCCTTCCGCAGTGTCATCTGCCACACGGAAGGCACGGCCACGGGCCAGGGG : 25864

   184 : GTCTATGAT  &lt;&lt;&lt;&lt; Intron 1 &lt;&lt;&lt;&lt;  CTTCTGCATGCCCAGCTGGCATGGCCCCA :   221
         |||||||||        404 bp        |||||||||||||||||||||||||||||
 25865 : GTCTATGATct..................acCTTCTGCATGCCCAGCTGGCATGGCCCCA : 26306

   222 : CGTAGAGT-GGNNTGGCGTCTCGGTGCTGGTCAGCGACACGTTGTCCTGGCTGGGCAGGT :   280
         |||||||| ||  |||||||||||||||||||||||||||||||||||||||||||||||
 26307 : CGTAGAGTGGGGGTGGCGTCTCGGTGCTGGTCAGCGACACGTTGTCCTGGCTGGGCAGGT : 26366

   281 : CCAGCTCCCGGAGGACCTGGGGCTTCAGCTTCCCGTAGCGCTGGCTGCAGTGACGGATGC :   340
         ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
 26367 : CCAGCTCCCGGAGGACCTGGGGCTTCAGCTTCCCGTAGCGCTGGCTGCAGTGACGGATGC : 26426

   341 : TCTTGCGCTGCCATTTCTGGGTGCTGTCACTGTCCTTGCTCACTCCAAACCAGTTCGGCG :   400
         ||||||||||||||||||||||||||||||||||||||||||||||||||||| ||||||
 26427 : TCTTGCGCTGCCATTTCTGGGTGCTGTCACTGTCCTTGCTCACTCCAAACCAG-TCGGCG : 26485

   401 : GTCCCC  &lt;&lt;&lt;&lt; Intron 2 &lt;&lt;&lt;&lt;  CTGCGGATGGTCTGTGTTGATGGACGTTTGGG :   438
         ||||||        898 bp        |||||||||||||||| |||||||||| | ||
 26486 : GTCCCCct..................acCTGCGGATGGTCTGTG-TGATGGACGTCT-GG : 27419

   439 : CTTTGCAGCACCGGCCGCC-GAGTTCATGG :   468
         | ||||||||||||||||| ||| ||||||
 27420 : CGTTGCAGCACCGGCCGCCGGAGCTCATGG : 27450



</pre>

<p> 

<hr>

<p> 

In the Ensembl CIGAR format the numbers and letters are switched,

and there are no gaps in the string.  So the above example in Ensembl would appear in a feature table in three rows with these CIGAR strings:

<p> 

<b>

<pre>
13M1I35M1I4M1I13M1D4M1I115M
37M1D164M1I12M
16M1I12M1I21M1D10M
</pre>

</b>

<p> 

<hr>

<strong>Related pages:</strong> <a href="/web/20050210133134/http://www.ensembl.org:80/Docs/wiki/html/EnsemblDocs/Exonerate.html">Exonerate</a>

<hr>


</td>
            </tr>
            <tr>
                <td><font size="2"><em>This page last edited on
                20 February 2003<br></em></font></td>
				<td align="right"><font size="2"><em><a href="/web/20050210133134/http://www.ensembl.org:80/cgi-bin/wiki/view.pl/EnsemblDocs/CigarFormat?rev=1.9">&lt;</a> Version:1.10 <br><a href="/web/20050210133134/http://www.ensembl.org:80/cgi-bin/wiki/view.pl/EnsemblDocs/CigarFormat?rev=1.10">printable version</a></em></font></td> 
            </tr>
        </table>
        <table border="0" cellpadding="0" cellspacing="0" width="100%" class="yellow1">
			<tr>
				<td>&nbsp;</td>
				<td width="1" bgcolor="#e2e2ff">&nbsp;</td>
				<td>&nbsp;</td>
			</tr>

            <tr>
                <td rowspan="2" width="50%">
				<ul>
                    <li><a href="/web/20050210133134/http://www.ensembl.org:80/Docs/wiki/html/EnsemblDocs/HomePage.html"><font size="2">EnsemblDocs Wiki HomePage</font></a><font size="2">
                        </font></li>
                    <li><a href="/web/20050210133134/http://www.ensembl.org:80/cgi-bin/wiki/search.pl/EnsemblDocs?text=Category&amp;options=topic"><font size="2">Wiki Categories</font></a></li>
                    <li><a href="/web/20050210133134/http://www.ensembl.org:80/cgi-bin/wiki/index.pl/EnsemblDocs"><font size="2">Index of EnsemblDocs Wiki</font></a></li>
                    <li><a href="/web/20050210133134/http://www.ensembl.org:80/cgi-bin/wiki/changes.pl/EnsemblDocs"><font size="2">Recent Changes to EnsemblDocs</font></a></li>
                    <li><a href="/web/20050210133134/http://www.ensembl.org:80/Docs/wiki/html/EnsemblDocs/FAQ.html"><font size="2">Wiki FAQ</font></a></li>
                </ul>
                </td>
                <td width="1" bgcolor="#e2e2ff">&nbsp;</td>
                <td valign="top"><ul>
                    <li><a href="/web/20050210133134/http://www.ensembl.org:80/cgi-bin/wiki/edit.pl/EnsemblDocs/CigarFormat"><font size="2">Edit this page</font></a></li>
                    <li><a href="/web/20050210133134/http://www.ensembl.org:80/cgi-bin/wiki/delete.pl/EnsemblDocs/CigarFormat"><font size="2">Delete this page</font></a></li>
                </ul>
                </td>
            </tr>
            <tr>
                <td width="1" bgcolor="#e2e2ff">&nbsp;</td>
                <td valign="top"><form action="/web/20050210133134/http://www.ensembl.org:80/cgi-bin/wiki/search.pl/EnsemblDocs/" method="POST">
                    <ul>
                        <li><font size="2">Search Wiki for:</font> <input type="text" size="20" name="text"></li>
                    </ul>
                </form>
                </td>
            </tr>
			<tr>
				<td>&nbsp;</td>
				<td width="1" bgcolor="#e2e2ff">&nbsp;</td>
				<td>&nbsp;</td>
			</tr>
        </table>
        </td>
    </tr>
</table>

<table border="0" cellpadding="0" cellspacing="0" width="100%" class="yellow1">
    <tr>
        <td bgcolor="#e2e2ff">&nbsp;</td>
        <td bgcolor="#e2e2ff">&nbsp;</td>
    </tr>
</table>


<!-- end rightbar -->
    </td>
  </tr>
</table>



<!-- begin footer -->

<br/>
<table border="0" width="100%%" cellpadding="0" cellspacing="0">
  <tr class="header2" valign="center">
    <td>
      &nbsp;<i>Date : Thu Feb 10 13:31:24 2005 </i>
    </td>
      <td align="center"><a href="https://web.archive.org/web/20050210133134/http://feb2005.archive.ensembl.org/Docs/wiki/html/EnsemblDocs/CigarFormat.html"><img valign="center" src="/web/20050210133134im_/http://www.ensembl.org:80/gfx/header/archive/jump2archive.gif" width="140" height="25" border="0"></a></td> 
    <td><img src="/web/20050210133134im_/http://www.ensembl.org:80/gfx/blank.gif" height="22" alt=""></td>
    <td align="RIGHT">
        <i><a href="/web/20050210133134/http://www.ensembl.org:80/helpdesk/">Help Desk / Suggestions</a></i>&nbsp;
    </td>
  </tr>
</table>
<br>
<!-- end footer -->

</body>
</html><!--
     FILE ARCHIVED ON 13:31:34 Feb 10, 2005 AND RETRIEVED FROM THE
     INTERNET ARCHIVE ON 15:23:36 Oct 29, 2018.
     JAVASCRIPT APPENDED BY WAYBACK MACHINE, COPYRIGHT INTERNET ARCHIVE.

     ALL OTHER CONTENT MAY ALSO BE PROTECTED BY COPYRIGHT (17 U.S.C.
     SECTION 108(a)(3)).
-->
