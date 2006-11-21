<?xml version="1.0"?>
<xsl:stylesheet version = "1.0" xmlns:xsl = "http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text"/>
<xsl:template match="/">
<xsl:text>Id dM Mass pepLen score expect z1 z2 z3 enzN enzC
</xsl:text>
<xsl:for-each select="//group[@type='model']">
<xsl:value-of select="concat(@id,'_',@z,' ')"/>
<xsl:value-of select="protein[1]/peptide/domain/@delta"/><xsl:text> </xsl:text>
<xsl:value-of select="@mh"/><xsl:text> </xsl:text>
<xsl:value-of select="string-length(protein[1]/peptide/domain/@seq)"/><xsl:text> </xsl:text>
<xsl:value-of select="protein[1]/peptide/domain/@hyperscore"/><xsl:text> </xsl:text>
<xsl:value-of select="protein[1]/peptide/domain/@expect"/><xsl:text> </xsl:text>
<xsl:if test="@z=1"><xsl:text>1 </xsl:text></xsl:if><xsl:if test="@z!=1"><xsl:text>0 </xsl:text></xsl:if>
<xsl:if test="@z=2"><xsl:text>1 </xsl:text></xsl:if><xsl:if test="@z!=2"><xsl:text>0 </xsl:text></xsl:if>
<xsl:if test="@z=3"><xsl:text>1 </xsl:text></xsl:if><xsl:if test="@z!=3"><xsl:text>0 </xsl:text></xsl:if>
<xsl:variable name="preC" select="substring(protein[1]/peptide/domain/@pre,string-length(protein[1]/peptide/domain/@pre),1)"/>
<xsl:variable name="pepN" select="substring(protein[1]/peptide/domain/@seq,1,1)"/>
<xsl:variable name="pepC" select="substring(protein[1]/peptide/domain/@seq,string-length(protein[1]/peptide/domain/@seq),1)"/>
<xsl:variable name="postN" select="substring(protein[1]/peptide/domain/@post,1,1)"/>
<xsl:variable name="trypN" select="$preC='[' or (($preC='K' or $preC='R') and $pepN!='P')"/>
<xsl:variable name="trypC" select="$postN=']' or (($pepC='K' or $pepC='R') and $postN!='P')"/>
<xsl:if test="$trypN"><xsl:text>1 </xsl:text></xsl:if><xsl:if test="$trypN!=true()"><xsl:text>0 </xsl:text></xsl:if>
<xsl:if test="$trypC"><xsl:text>1 </xsl:text></xsl:if><xsl:if test="$trypC!=true()"><xsl:text>0 </xsl:text></xsl:if>
<xsl:text>
</xsl:text>
</xsl:for-each>
</xsl:template></xsl:stylesheet>

