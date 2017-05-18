Attribute VB_Name = "NewMacros"
Sub wrap_table()
Attribute wrap_table.VB_ProcData.VB_Invoke_Func = "Normal.NewMacros.Macro3"
'
' wrap_table Macro
'
'
    Dim t As Table
    For Each t In ActiveDocument.Tables
        t.Rows.WrapAroundText = True
        t.Rows.HorizontalPosition = wdTableLeft
        t.Rows.AllowOverlap = False
    Next
    Selection.ParagraphFormat.Alignment = wdAlignParagraphJustify
End Sub
