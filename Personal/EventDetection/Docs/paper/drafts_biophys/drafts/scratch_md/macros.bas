Attribute VB_Name = "NewMacros"
Sub wrap_table()
Attribute wrap_table.VB_ProcData.VB_Invoke_Func = "Normal.NewMacros.Macro3"
'
' wrap_table Macro
'
'
    ' Make an index array
    '
    Dim i As Integer
    i = 0
    Dim myarray(1 To 3) As Integer
    Dim t As Table
    For Each t In ActiveDocument.Tables
        i = i + 1
        t.Rows.WrapAroundText = True
        t.Rows.HorizontalPosition = wdTableLeft
        t.Rows.VerticalPosition = InchesToPoints(1)
        t.Rows.RelativeVerticalPosition = wdRelativeVerticalPositionPage
        t.Rows.AllowOverlap = False
    Next
    Selection.ParagraphFormat.Alignment = wdAlignParagraphJustify
End Sub
